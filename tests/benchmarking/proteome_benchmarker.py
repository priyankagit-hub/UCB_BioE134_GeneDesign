import os
import traceback
import csv
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.seq_utils.translate import Translate

def parse_fasta(fasta_file):
    """
    Parses the FASTA file to extract gene names and protein sequences.
    
    Parameters:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        dict: A dictionary where keys are gene names and values are protein sequences.
    """
    sequences = {}
    current_gene = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # If we encounter a new gene, save the previous one
                if current_gene:
                    sequences[current_gene] = ''.join(current_sequence)
                # Reset for the new gene
                current_gene = line.split('|')[2].split(' ')[0]  # Extract gene name
                current_sequence = []
            else:
                # Append the sequence lines
                current_sequence.append(line)

        # Save the last sequence
        if current_gene:
            sequences[current_gene] = ''.join(current_sequence)
    
    return sequences

def benchmark_proteome(fasta_file):
    """
    Benchmarks the proteome using TranscriptDesigner.
    
    Parameters:
        fasta_file (str): Path to the FASTA file.
    
    Logs:
        The gene name, protein sequence, and result (or error).
    """
    # Initialize TranscriptDesigner
    designer = TranscriptDesigner()
    designer.initiate()

    # Parse the FASTA file
    proteome = parse_fasta(fasta_file)

    # Prepare to log results
    successful_results = []
    error_results = []

    for gene, protein in proteome.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {protein[:30]}...")  # Show first 30 aa for brevity
            
            # Run the TranscriptDesigner on the protein sequence
            ignores = set()
            transcript = designer.run(protein, ignores)
            
            # Append successful result
            successful_results.append({
                'gene': gene,
                'protein': protein,
                'transcript': transcript
            })
        
        except Exception as e:
            # Capture error and stack trace
            error_results.append({
                'gene': gene,
                'protein': protein,
                'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"
            })
    
    return successful_results, error_results

def analyze_errors(error_results):
    """
    Analyze the errors and write them to a text file.
    
    Parameters:
        error_results (list): List of error dictionaries.
    """
    error_summary = {}
    
    # Collect statistics about the types of errors
    for error in error_results:
        error_message = error['error'].split("\n")[0]
        error_summary[error_message] = error_summary.get(error_message, 0) + 1

    # Write the error analysis to a text file
    with open('error_summary.txt', 'w') as f:
        for error_message, count in error_summary.items():
            f.write(f"{error_message}: {count} occurrences\n")

        f.write("\nDetailed errors:\n")
        for error in error_results:
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")

def validate_transcripts(successful_results):
    """
    Validate the successful transcripts using various checkers.
    
    Parameters:
        successful_results (list): List of successful result dictionaries.
    
    Logs:
        Validation results, including any failures.
    """
    # Initialize ForbiddenSequenceChecker, PromoterChecker, and Translate
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()

    promoter_checker = PromoterChecker()
    promoter_checker.initiate()

    translator = Translate()
    translator.initiate()

    validation_failures = []

    for result in successful_results:
        # Extract CDS (Coding DNA Sequence)
        cds = ''.join(result['transcript'].codons)

        try:
            # Check CDS Translation and Completeness
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")

            original_protein = result['protein']
            translated_protein = translator.run(cds)

            if original_protein != translated_protein:
                raise ValueError(f"Translation mismatch: Original {original_protein}, Translated {translated_protein}")

            # Check CDS for Start and Stop codons
            if not (cds.startswith(("ATG", "GTG", "TTG")) and cds.endswith(("TAA", "TGA", "TAG"))):
                raise ValueError("CDS does not start with a valid start codon or end with a valid stop codon.")

        except ValueError as e:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Translation or completeness error: {str(e)}"
            })
            continue  # Skip further checks for this transcript due to translation issues

        # Rebuild full transcript DNA sequence (RBS + CDS)
        transcript_dna = result['transcript'].rbs.utr.upper() + cds

        # Check for Minimal Hairpins
        passed_hairpin, hairpin_string = hairpin_checker(transcript_dna)
        if not passed_hairpin:
            # Format the hairpin string for TSV output, ensuring there are no unnecessary quotes
            formatted_hairpin = hairpin_string.replace('\n', ' ').replace('"', "'")
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Hairpin detected: {formatted_hairpin}"
            })

        # Check for Forbidden Sequences
        passed_forbidden, forbidden_site = forbidden_checker.run(transcript_dna)
        if not passed_forbidden:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Forbidden sequence: {forbidden_site}"
            })

        # Check for Internal Promoters
        passed_promoter, found_promoter = promoter_checker.run(transcript_dna)
        if not passed_promoter:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Constitutive promoter detected: {found_promoter}" if found_promoter else "Constitutive promoter detected"
            })


    # Output the validation failures as TSV
    with open('validation_failures.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene', 'protein', 'cds', 'site'])
        for failure in validation_failures:
            writer.writerow([failure['gene'], failure['protein'], failure['cds'], failure['site']])

def run_benchmark(fasta_file):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, and validating.
    
    Parameters:
        fasta_file (str): Path to the FASTA file.
    """
    # Benchmark the proteome
    successful_results, error_results = benchmark_proteome(fasta_file)

    # Analyze and log errors
    if error_results:
        print(f"Analyzing {len(error_results)} errors...")
        analyze_errors(error_results)

    # Validate the successful transcripts
    if successful_results:
        print(f"Validating {len(successful_results)} successful transcripts...")
        validate_transcripts(successful_results)

if __name__ == "__main__":
    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    
    # Run the complete benchmark
    run_benchmark(fasta_file)
