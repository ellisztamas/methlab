import pandas as pd
import numpy as np

import methlab as ml

def test_estimate_beta_parameters():
    """Test that function recovers inputs.
    """
    a,b = 17, 22
    mu = a / (a+b)
    sigma = (a*b) / ( (a+b)**2 * (a+b+1) )
    
    hat = ml.estimate_beta_parameters(mu, sigma)
    assert hat == (a,b)

read_counts = pd.read_csv("tests/test_data/gene_read_counts.csv")

def test_methylation_state():
    # Check that when error rates are very low, almost everything looks TE-like methylated
    # There is one feature where methylated is *really* zero
    ab_errors = (0.1, 220)
    low_errors = ml.methylation_state(read_counts, ab_errors, hard_calls = True)
    assert np.sum(low_errors['call'] == "TE-like") == 9

    # Check that when error rates are very low, almost everything looks TE-like methylated
    # There is one feature where methylated is *really* zero
    ab_errors = (220,17)
    high_errors = ml.methylation_state(read_counts, ab_errors, hard_calls=True)
    assert np.sum(high_errors['call'] == "unmethylated") == 9