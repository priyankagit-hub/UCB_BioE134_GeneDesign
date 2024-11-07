import pytest
from genedesign.checkers.gc_checker import GCContentChecker

@pytest.fixture
def gc_content_checker():
    """
    Fixture to initialize the GCContentChecker with a specific GC range.
    """
    return GCContentChecker(min_gc=0.45, max_gc=0.55)

def test_gc_within_range(gc_content_checker):
    """
    Test case for a sequence within the GC range, which should pass.
    """
    sequence = "GCGCATGCGC"  # 50% GC content
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == True
    assert 0.45 <= gc_content <= 0.55

def test_gc_below_min_range(gc_content_checker):
    """
    Test case for a sequence below the GC range, which should fail.
    """
    sequence = "ATATATATAT"  # Low GC content (0% GC)
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == False
    assert gc_content < 0.45

def test_gc_above_max_range(gc_content_checker):
    """
    Test case for a sequence above the GC range, which should fail.
    """
    sequence = "GCGCGCGCGC"  # High GC content (100% GC)
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == False
    assert gc_content > 0.55

def test_gc_at_lower_bound(gc_content_checker):
    """
    Edge case test for a sequence at the lower bound of GC range, which should pass.
    """
    sequence = "ATGCATGC"  # 50% GC content
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == True
    assert 0.45 <= gc_content <= 0.55

def test_gc_at_upper_bound(gc_content_checker):
    """
    Edge case test for a sequence just at the upper bound of the GC range, which should pass.
    """
    sequence = "GCGCATGC"  # 62.5% GC content, fails with this threshold
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == False
    assert gc_content > 0.55

def test_empty_sequence(gc_content_checker):
    """
    Test case for an empty sequence, which should fail by default with 0% GC content.
    """
    sequence = ""
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == False
    assert gc_content == 0.0

def test_mixed_gc_content(gc_content_checker):
    """
    Test case for a mixed GC content sequence, close to mid-range.
    """
    sequence = "ATGCGCGAATCG"  # Around 50% GC content
    passed, gc_content = gc_content_checker.run(sequence)
    
    assert passed == True
    assert 0.45 <= gc_content <= 0.55
