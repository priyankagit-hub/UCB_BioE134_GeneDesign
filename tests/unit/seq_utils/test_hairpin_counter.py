
import pytest
from genedesign.seq_utils.hairpin_counter import hairpin_counter

def test_hairpin_counter():
    # Test cases for different sequences
    assert hairpin_counter("AAAAAAAAAAAAAAAAAAAAAAAAAAA") == 0, "Test case 1 failed: No hairpins expected"
    assert hairpin_counter("AAAAACCCAAAAAAAAAAGGGAAAAAA") == 0, "Test case 2 failed: No hairpin expected"
    assert hairpin_counter("AAAAACCCAAAAAAAAAGGGAAAAAA") == 1, "Test case 3 failed"
    assert hairpin_counter("AAAAACCCCAAAAAAAAGGGGAAAAAA") == 3, "Test case 4 failed"
    assert hairpin_counter("AAAAACACGAAAAAAAACGTGAAAAAA") == 1, "Test case 5 failed"
    assert hairpin_counter("AAAACCCCCAAAAAAAAGGGGGAAA") == 3, "Test case 6 failed"
    assert hairpin_counter("AAAACCCCCAAAAAAAGGGGGAAA") == 6, "Test case 7 failed"

if __name__ == "__main__":
    pytest.main()