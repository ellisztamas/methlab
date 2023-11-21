import pandas as pd
import numpy as np
from scipy.stats import beta, binom
from methlab._alogsumexp import alogsumexp

def methylation_state(read_counts:pd.DataFrame, ab_errors:tuple, return_probabilities:bool=False, hard_calls:bool=False):
    """
    Likelihoods that windows or features of a genome are in one several distinct
    states.

    This takes a table of observed counts of converted and unconverted reads 
    within windows or features in each sequence context and calculates
    likelihoods that each is in one of several distinct methylation 'states'.
    The current implementation assesses the evidence for three states:
    - Unmethylated (no methylated in any sequence context).
    - Gene-body-like methylation (CG, but no CHG or CHH methylation).
    - TE-like methylation (methylation in all three contexts).

    This compares the evidence that observed methylated read counts patterns were
    generated from:
    - the expected distribution of non-conversion errors, modelled as a
        beta distribution
    - an additional process, modelled as the cumulative distribution of the
        same beta distribution.

    Parameters
    ==========
    read_counts: dataframe
        Pandas dataframe with (at least) columns 'id' giving a name for each 
        feature,  'context' giving sequence context (CG, CHG, CHH, total),
        'unconverted' and 'converted' (integer counts of unconverted and
        converted reads mapping in each context; this can be generated from
        methylation_over_features() or methylation_in_windows() from
        CytosineCoverageFile).
    ab_errors: tuple
        Tuple of length 2 giving the a and b shape parameters of the beta
        distribution describing variation in non-conversion errors.
    return_probabilities: bool
        If True, likelihoods for each state are converted to probabilities that
        sum to one.
    hard_calls: bool
        If True, returns a column giving the most-likely state. Defaults to
        False.
    
    Returns
    =======
    A dataframe giving ID, coverage, and (log) likelihoods that each
    region/window is unmethylated, CG-methylated only, or TE-like methylated. If
    `return_probabilities` is True, likelihoods are reported as probabilities
    summing to one.
    """

    # Mean methylation
    read_counts['coverage'] = read_counts['unconverted'] + read_counts['converted']
    read_counts['theta']    = read_counts['unconverted'] / read_counts['coverage']

    read_counts.loc[read_counts['theta'] == 0, 'theta'] = 1e-12
    read_counts.loc[read_counts['theta'] == 1, 'theta'] = 1 - 1e-12


    # Log binomial probabilities of the data, given the observed means
    log_prob_binom = binom.logpmf(
            k = read_counts['unconverted'],
            n = read_counts['coverage'],
            p = read_counts['theta']
            )
    # Log probabilities of the observed mean under the error distribution, or a distribution greater than that.
    prob_means = {
        'errors' : log_prob_binom + beta.logpdf(
            x = read_counts['theta'],
            a = ab_errors[0],
            b = ab_errors[1]
            ),
        'alternative' : log_prob_binom + beta.logcdf(
            x = read_counts['theta'],
            a = ab_errors[0],
            b = ab_errors[1]
            )
    }

    # To simplify syntax later, make a dictionary of index positions for each context
    ix = {
        'CG'  : read_counts['context'] == "CG",
        'CHG' : read_counts['context'] == "CHG",
        'CHH' : read_counts['context'] == "CHH",
    }

    # Likelihoods of generating the data if:
        # 1. that methylation is all due to non-conversion errors
        # 2. non-CG methylation is due to non-conversion errors, but CG is not
        # 3. methylation in all contexts is greater than expected from non-conversion errors
    state_likelihoods = {
        'unmethylated' : prob_means['errors'][ix['CG']]      + prob_means['errors'][ix['CHG']]      + prob_means['errors'][ix['CHH']],
        'CG-only'      : prob_means['alternative'][ix['CG']] + prob_means['errors'][ix['CHG']]      + prob_means['errors'][ix['CHH']],
        'TE-like'      : prob_means['alternative'][ix['CG']] + prob_means['alternative'][ix['CHG']] + prob_means['alternative'][ix['CHH']]
    }

    # Normalise likelihoods to sum to one.
    if return_probabilities:
        row_sums = alogsumexp(
            np.array([
                state_likelihoods['unmethylated'],
                state_likelihoods['CG-only'],
                state_likelihoods['TE-like']
            ]),
            axis=0
        )
        state_likelihoods = { k : np.exp(v - row_sums) for k,v in state_likelihoods.items() }
    
    # Output file with only IDs and coverage
    output = read_counts.loc[read_counts['context'] == "total"][['id', 'ncytosines', 'coverage', 'theta']].reset_index()
    output = pd.concat([output, pd.DataFrame(state_likelihoods)], axis =1)
    output = output.drop(['index'], axis=1)

    # Make a call about methylation status
    if hard_calls:
        a = np.array([
            state_likelihoods['unmethylated'],
            state_likelihoods['CG-only'],
            state_likelihoods['TE-like']
        ]).T
        max_ix = np.argmax(a, axis =1) + 1
        max_ix[np.isnan(a).sum(axis=1) == 3] = 0

        output['call'] = [ ['NA', 'unmethylated', 'CG-only', "TE-like"][i] for i in max_ix ]

    return(output)
