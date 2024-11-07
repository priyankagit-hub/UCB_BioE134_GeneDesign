class GCContentChecker:
    def __init__(self, min_gc=0.45, max_gc=0.55):
        """
        Initialize with minimum and maximum GC content thresholds.
        """
        self.min_gc = min_gc
        self.max_gc = max_gc

    def run(self, sequence: str) -> tuple[bool, float]:
        """
        Check if the GC content of the sequence is within the defined bounds.

        Parameters:
            sequence (str): The DNA sequence to check.

        Returns:
            tuple: (bool, float) where the boolean indicates if the sequence passes
            and the float is the calculated GC content.
        """
        if not sequence:
            return False, 0.0  # Handle empty sequences

        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence)

        is_within_bounds = self.min_gc <= gc_content <= self.max_gc
        return is_within_bounds, gc_content

