import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.gc_checker import GCContentChecker

class TranscriptDesigner:
    def __init__(self, codon_usage_file="genedesign/data/codon_usage.txt", seed=42):
        random.seed(seed)
        
        # Initialize components
        self.rbs_chooser = RBSChooser()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        self.codon_checker = CodonChecker()
        self.gc_checker = GCContentChecker()
        
        # Parameters
        self.codon_usage_file = codon_usage_file
        self.window_size = 3
        self.downstream_size = 8
        self.max_attempts = 50
        self.preamble_codons_count = 8

        # Codon weights and precomputed lists
        self.codon_weights = self.load_codon_usage(self.codon_usage_file)
        self.weighted_codon_lists = self.precompute_weighted_codon_lists()

    def initiate(self):
        """Initialize all checkers."""
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()
        self.codon_checker.initiate()
        self.rbs_chooser.initiate()

    def load_codon_usage(self, filepath):
        """Load codon usage frequencies."""
        codon_weights = {}
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    codon, aa, freq = parts[0], parts[1], float(parts[2])
                    if aa not in codon_weights:
                        codon_weights[aa] = []
                    codon_weights[aa].append((codon, freq))
        return codon_weights

    def precompute_weighted_codon_lists(self):
        """Precompute weighted codon lists for each amino acid."""
        weighted_codon_lists = {}
        for aa, codons in self.codon_weights.items():
            weighted_list = []
            for codon, freq in codons:
                count = max(int(round(freq * 100)), 1)  
                weighted_list.extend([codon] * count)
            weighted_codon_lists[aa] = weighted_list
        return weighted_codon_lists

    def select_random_codon(self, aa):
        """Select a random codon for the given amino acid."""
        return random.choice(self.weighted_codon_lists[aa])

    def segment_passes_all_checks(self, segment, codons):
        """Check if a segment passes all checks without scoring penalties."""
        # Forbidden Sequence Checker
        forbidden_passed, _ = self.forbidden_checker.run(segment)
        if not forbidden_passed:
            return False

        # Internal Promoter Checker
        promoter_passed, _ = self.promoter_checker.run(segment)
        if not promoter_passed:
            return False

        # Hairpin Checker
        hairpin_passed, _ = hairpin_checker(segment)
        if not hairpin_passed:
            return False

        # Codon Usage Checker
        codons_above_board, _, _, _ = self.codon_checker.run(codons)
        if not codons_above_board:
            return False

        # GC Content Checker
        gc_passed, _ = self.gc_checker.run(segment)
        if not gc_passed:
            return False

        return True

    def score_segment(self, segment, codons):
        """Score a DNA sequence segment based on various checkers."""
        score = 0
        forbidden_passed, _ = self.forbidden_checker.run(segment)
        score += 30 if forbidden_passed else -50

        promoter_passed, _ = self.promoter_checker.run(segment)
        score += 30 if promoter_passed else -50

        hairpin_passed, _ = hairpin_checker(segment)
        score += 50 if hairpin_passed else -100

        codons_above_board, diversity, rare_codons, cai = self.codon_checker.run(codons)
        if codons_above_board:
            score += int(diversity * 50) + int(cai * 100) - rare_codons * 10
        else:
            score -= 100

        gc_passed, _ = self.gc_checker.run(segment)
        score += 20 if gc_passed else -30

        return score

    def max_score(self):
        """Calculate the maximum possible score based on scoring weights."""
        return 50 + 30 + 40 + 50 + 100 + 20

    def monte_carlo_window(self, window_peptide, codons_so_far, downstream_peptide):
        """Finds the best codon sequence for a window."""
        preamble_seq = ''.join(codons_so_far[-self.preamble_codons_count:] if codons_so_far else [])

        best_codons, best_score = None, -float('inf')

        for _ in range(self.max_attempts):
            window_codons = [self.select_random_codon(aa) for aa in window_peptide]
            downstream_codons = [self.select_random_codon(aa) for aa in downstream_peptide[:self.downstream_size]]
            full_seq = preamble_seq + ''.join(window_codons) + ''.join(downstream_codons)
            segment_codons = codons_so_far[-self.preamble_codons_count:] + window_codons

            # Check if this segment passes all criteria
            if self.segment_passes_all_checks(full_seq, segment_codons):
                return window_codons  # Immediately accept if it passes all checks

            # Otherwise, score the segment
            score = self.score_segment(full_seq, segment_codons)
            if score > best_score:
                best_score = score
                best_codons = window_codons

        # Return the highest scoring option if none fully passed
        if best_codons is None:
            raise RuntimeError("Unable to find valid codon sequence.")
        
        return best_codons

    def run(self, peptide, ignores=set()):
        """Designs a transcript for the peptide sequence."""
        if self.codon_weights is None:
            raise RuntimeError("TranscriptDesigner not initiated. Please call 'initiate()' before 'run()'.")

        codons = []
        current_index = 0

        while current_index < len(peptide):
            # Define the window peptide
            window_peptide = peptide[current_index:current_index + self.window_size]
            downstream_peptide = peptide[current_index + self.window_size:current_index + self.window_size + self.downstream_size]

            # Generate codons for the current window
            window_codons = self.monte_carlo_window(window_peptide, codons, downstream_peptide)
            codons.extend(window_codons)

            # Slide the window forward by 3 codons
            current_index += self.window_size

        # Handle any remaining amino acids not covered by the window
        remaining_aas = peptide[current_index:]
        for aa in remaining_aas:
            codons.append(self.select_random_codon(aa))

        # Add stop codon
        codons.append('TAA')  # You can choose the most frequent stop codon if preferred

        # Create the complete CDS
        cds = ''.join(codons)

        # Select RBS (not using 'ignores' in this simplified version)
        selected_rbs = self.rbs_chooser.run(cds, ignores)

        # Create the Transcript object
        return Transcript(selected_rbs, peptide, codons)