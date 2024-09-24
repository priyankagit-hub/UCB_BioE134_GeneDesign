import pytest
from genedesign.checkers.hairpin_counter import hairpin_counter

def test_no_hairpin():
    sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAA"
    result = hairpin_counter(sequence)
    assert result == 0, "Expected no hairpins, but got a non-zero count."

def test_hairpin_detected():
    sequence = "CCCCCTTTCCCCCCAAACCCCCC"  # Example of a weak hairpin
    result = hairpin_counter(sequence)
    assert result > 0, "Expected a hairpin, but got zero."

def test_strong_hairpin():
    sequence = "AAAAACCCCCAAAAAAAAGGGGGAAA"  # Example of a strong hairpin
    weak_sequence = "AAAAACCCAAAAAAAAAGGGAAAAAA"  # Example of a weaker hairpin
    
    strong_result = hairpin_counter(sequence)
    weak_result = hairpin_counter(weak_sequence)
    
    assert strong_result > weak_result, "Expected stronger hairpin to have a higher count than weaker one."
    assert strong_result > 0, "Expected a hairpin in the strong sequence, but got zero."
    assert weak_result > 0, "Expected a hairpin in the weak sequence, but got zero."

@pytest.mark.parametrize("sequence", [
    "GGCTAATTTAGCCATTAAGGCTAATAGGCTAA",  # Example: may or may not form a hairpin
    "AAAAACACGAAAAAAAACGTGAAAAAA",       # Potential hairpin
    "AAAACCCCCAAAAAAAAGGGGGAAA",         # Strong hairpin
])
def test_varied_hairpins(sequence):
    result = hairpin_counter(sequence)
    assert result >= 0, "Hairpin count should never be negative."

