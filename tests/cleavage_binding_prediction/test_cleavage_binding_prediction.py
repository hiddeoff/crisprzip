import pytest
import pandas as pd
from crisprzip.kinetics import *

PRECISION = 5E-4

# Load expected output
DATA_AVG = pd.read_csv("tests/cleavage_binding_prediction/expected_output_avg.csv",
                     dtype={'mm_pattern': str})
DATA_SEQ = pd.read_csv("tests/cleavage_binding_prediction/expected_output_seq.csv")


@pytest.fixture
def average_searcher() -> Searcher:
    """Loads average-based parameter set."""
    return load_landscape('average_params')


@pytest.fixture
def sequence_searcher() -> BareSearcher:
    """Loads sequence-based parameter set."""
    return load_landscape('sequence_params')


# noinspection PyTypeChecker
@pytest.mark.parametrize("mm_pattern, exptype, k_on, time, fraction", DATA_AVG.values)
def test_average_model(mm_pattern, exptype, k_on, time, fraction, average_searcher):
    mm_pattern = MismatchPattern.from_string(mm_pattern)
    complex_obj = average_searcher.probe_target(mm_pattern)
    if exptype == 'cleavage':
        f_clv = complex_obj.get_cleaved_fraction(time=time, on_rate=k_on)
        assert f_clv == pytest.approx(fraction, abs=PRECISION)
    elif exptype == 'binding':
        f_clv = complex_obj.get_bound_fraction(time=time, on_rate=k_on)
        assert f_clv == pytest.approx(fraction, abs=PRECISION)
    else:
        raise ValueError(f"Unknown exptype '{exptype}'")


# noinspection PyTypeChecker
@pytest.mark.parametrize("protospacer, targetseq, exptype, k_on, time, fraction", DATA_SEQ.values)
def test_sequence_model(protospacer, targetseq, exptype, k_on, time, fraction, sequence_searcher):
    complex_obj = sequence_searcher.probe_sequence(protospacer=protospacer, target_seq=targetseq)
    if exptype == 'cleavage':
        f_clv = complex_obj.get_cleaved_fraction(time=time, on_rate=k_on)
        assert f_clv == pytest.approx(fraction, abs=PRECISION)
    elif exptype == 'binding':
        f_clv = complex_obj.get_bound_fraction(time=time, on_rate=k_on)
        assert f_clv == pytest.approx(fraction, abs=PRECISION)
    else:
        raise ValueError(f"Unknown exptype '{exptype}'")
