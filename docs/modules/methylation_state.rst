===========================
Discrete methylation states
===========================


We wish to infer whether some region or window of DNA is in one of several discrete states.
For example, one model of plant methylation is that DNA is

- Unmethylated (no converted reads in any sequence context)
- Gene-body-like, or CG-only methylated (converted reads on CG sites, but not CHG or CHH sites)
- TE-like methylated (converted sites in CG, CHG and CHH contexts)

This is challenging when there are substantial non-conversion errors, especially when these vary across the genome.
Fortunately, for this task we do not need to accurately estimate true methylation level.
Instead, we need only ask whether the data are consistent with non-conversion errors alone, or with the cumulative effect of errors plus some additional process.

Note that it is also possible that there may be false conversion errors (cytosines that are observed as converted to thymine, but are truly unmethylated).
The current implementation can handle non-conversion errors only.

Data
====

Functions from `CytosineCoverageFile` will give you data that looks something like this, summarising feature name, sequence context, numbers of converted and unconverted reads, and total number of cytosines.

.. parsed-literal::
            id context  unconverted  converted  ncytosines
    0   feature0      CG     2            14          20
    1   feature0     CHG     0             4           8
    2   feature0     CHH     6            78         106
    3   feature0   total     8            96         134
    4   feature1      CG     2            14          20
    5   feature1     CHG     0             4           8
    6   feature1     CHH     6            78         106
    7   feature1   total     8            96         134
    8   feature2      CG     2            13          18
    9   feature2     CHG     0             4           6
    10  feature2     CHH     3            53          64
    11  feature2   total     5            70          88

Errors as a beta distribution
=============================

We assume non-conversion errors vary across the genome, such that we can't know the precise error rate at any single region.
However, we do have a good idea of the distribution of error rates across the genome, and model this as a beta distribution.
Whereas a binomial distribution describes the probability of observing binary data given some mean, a beta distribution describes the uncertainty in that mean.

A beta distribution is described by shape parameters *a* and *b*.
These can be estimated by method-of-moments based on the mean and variance of a samples from a beta distribution.
This code imports an example file from the output of `CytosineCoverageFile.methylation_over_features()` estimated in windows from control DNA known to be unmethylated, and estimates mean non-conversion in each window.

.. code-block:: python

    import pandas as pd
    import methlab as ml
    read_counts = pd.read_csv("tests/test_data/gene_read_counts.csv")
    # Subset to include only rows showing total counts for all cytosines
    read_counts = read_counts.loc[read_counts['context'] == "total"]
    # Estimate the proportion of unconverted reads
    read_counts['n'] = read_counts['unconverted'] + read_counts['converted']
    read_counts['theta'] = read_counts['unconverted'] / read_counts['n']

Estimate the mean error rates within each window, and the variance between windows

.. code-block:: python

    mu = read_counts['theta'].mean()
    sigma = read_counts['theta'].var()

This ignores differences in coverage and number of cytosines between windows.
You might want to include those as weights.

``estimate_beta_parameters()`` uses these values to estimate the shape parameters of the beta distribution.

.. code-block:: python
    
    ab_errors = ml.estimate_beta_parameters(mu, sigma)

It returns a tuple of two numbers corresponding to the *a* and *b* parameters.

Methylation state
=================

We can use the table of read counts and the shape of the beta distribution to calculate likelihoods that each feature/window in the table is unmethylated, gene-body-like methylated or TE-like methylated.

.. code-block:: python

    ml.methylation_state(read_counts, ab_errors)

Here is the output as a Pandas dataframe:

.. parsed-literal:: 

              id  coverage  unmethylated    CG-only    TE-like
    0  AT1G01010       746     -5.562562  -5.935396  -8.910764
    1  AT1G01020       840     -4.486843  -5.526963  -9.010684
    2  AT1G01030       882     -4.962632  -6.173149  -9.211579
    3  AT1G01040      2068     -8.321224  -7.392753 -10.439014
    4  AT1G01046        25     46.701561  19.947363 -33.561034
    5  AT1G01050       634     -3.166360  -3.806465  -8.671619
    6  AT1G01060      1156     -8.514358  -8.252310  -9.607966
    7  AT1G01070      1052     -6.517405  -7.524795  -9.311916
    8  AT1G01073        34    -29.364506  -5.335205  -3.258661
    9  AT1G01080       643     -5.348341  -6.894794  -8.887364

The last three columns give (log) likelihoods that each gene is in each of the three states.
If you include the argument ``return_probabilities=True`` log likelihoods are converted to probabilities that sum to one.
If you include the argument ``hard_calls=True`` an additional column is added giving the most likely state.
This can be convenient, but be aware that generally you probably want to *quantify* the evidence for each hypothesis rather than make hard calls like this, which is why it defaults to ``False``.

