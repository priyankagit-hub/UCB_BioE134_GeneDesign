import csv
import pandas as pd  # type: ignore
from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.Translate import Translate
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from typing import Set
import logging

class RBSChooser:
    """
    A class to select the best Ribosome Binding Site (RBS) for a given coding sequence (CDS),
    using a default merged data file unless specified otherwise.
    """
    rbs_options = set()  # Class variable to hold all RBS options

    def __init__(self, merged_data_file: str = "genedesign/data/merged_data.csv"):
        self.translator = Translate()
        self.translator.initiate()
        self.merged_data_file = merged_data_file

    def initiate(self):
        """
        Populates the RBS options from the merged data CSV file, treating the first column as locus_tag.
        """
        df = pd.read_csv(self.merged_data_file, index_col=0)
        
        # Iterate over rows and create RBSOptions
        for locus_tag, row in df.iterrows():
            utr = row['UTR']
            cds = row['CDS']
            gene_name = row['gene']
            
            # Translate the first 18 bases (first 6 amino acids)
            first_six_aas = self.translator.run(cds[:18])
            rbs_option = RBSOption(
                utr=utr,
                cds=cds,
                gene_name=gene_name,
                first_six_aas=first_six_aas
            )
            RBSChooser.rbs_options.add(rbs_option)

    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS using available RBS options.
        
        Parameters:
            cds (str): The coding DNA sequence for which an RBS is to be selected.
            ignores (Set[RBSOption]): A set of RBS options to ignore during selection.

        Returns:
            RBSOption: The selected RBSOption object that best fits the given CDS.
        """
        # Exclude ignored RBS options
        valid_rbs_options = [rbs for rbs in RBSChooser.rbs_options if rbs not in ignores]

        # Evaluate secondary structure (hairpin count) for each valid option
        rbs_hairpin_scores = []
        for rbs in valid_rbs_options:
            combined_seq = rbs.utr + cds[:30]
            hairpin_count, _ = hairpin_counter(combined_seq)
            rbs_hairpin_scores.append((rbs, hairpin_count))

        # Peptide similarity comparison (edit distance between first six amino acids)
        rbs_scores = []
        first_six_aas_cds = self.translator.run(cds[:18])
        for rbs, hairpin_count in rbs_hairpin_scores:
            edit_distance = calculate_edit_distance(rbs.first_six_aas, first_six_aas_cds)
            # Final score combines hairpin count and edit distance
            final_score = hairpin_count + edit_distance
            rbs_scores.append((rbs, final_score))

        # Sort RBS options by final score (lower score is better)
        rbs_scores.sort(key=lambda x: x[1])

        if not rbs_scores:
            raise ValueError("No valid RBS options found after scoring.")
        
        # Return the RBS with the best (lowest) score

        return rbs_scores[0][0]  

