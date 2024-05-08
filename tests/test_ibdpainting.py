import numpy as np
import methlab as ml
import pytest
from warnings import warn

from methlab.ibdpainting import import_vcfdistance

input = 'tests/test_data/panel_to_test.vcf.gz'
reference = 'tests/test_data/reference_panel.vcf.gz'
chr1 = 'tests/test_data/reference_panel_chr1.vcf.gz'

def test_import_vcfdistance_gives_right_output():
    """
    Check that import_vcfdistance gives the right answers when it should.
    """
    out = import_vcfdistance(
        input = input,
        reference = reference,
        sample_name = '1158_2'
    )
    real_names = np.array(['1158_2', '1158', '6024', '6184', '8249'])
    assert all(
        [ x == y for x,y in zip(out.samples, real_names) ]
    )

    assert len(out.chr) == 551
    assert len(out.pos) == 551
    assert out.geno.shape == (551, 5, 2)

def test_import_vcfdistance_fails_if_missing_sample():
    """
    Check that import_vcfdistance fails if the sample name is not in the input VCF
    """
    with pytest.raises(Exception):
        import_vcfdistance(
            input = input,
            reference = reference,
            sample_name = 'not_a_real_sample_name'
        )

def test_import_vcfdistance_fails_if_contigs_dont_match():
    with pytest.raises(Exception):
        import_vcfdistance(
            input = input,
            reference = chr1,
            sample_name = '1158_2'
        )

def test_split_into_windows_functions():
    vcfd = import_vcfdistance(
            input = input,
            reference = reference,
            sample_name = '1158_2'
        )
    split_vcfd = vcfd.split_into_windows(1000)
    assert all( split_vcfd['Chr1:0-1000'].pos >= 0 )
    assert all( split_vcfd['Chr1:0-1000'].pos < 1000 )
    assert all(split_vcfd['Chr1:0-1000'].chr == "Chr1")
    assert len(split_vcfd['Chr1:0-1000'].geno.shape) == 3
    # Check you get only one window per chr if window size >> chr length
    assert len(vcfd.split_into_windows(1000000)) == 2

def test_pairwise_distance_works():
    """
    There are four accessions in the reference VCF.
    Test each against the whole panel, and check that one of them comes out as
    identical in each case.
    """
    # 1158
    check_1158 = import_vcfdistance(input = reference, reference = reference,
                sample_name= '1158'
        ).pairwise_distance()

    assert check_1158[0] == 0
    assert all(check_1158[1:] > 0)

    # 6024
    check_6024 = import_vcfdistance(input = reference, reference = reference,
                sample_name= '6024'
        ).pairwise_distance()

    assert check_6024[1] == 0
    assert all(check_6024[[0,2,3]] > 0)

    # 6184
    check_6184 = import_vcfdistance(input = reference, reference = reference,
                sample_name= '6184'
        ).pairwise_distance()

    assert check_6184[2] == 0
    assert all(check_6184[[0,1,3]] > 0)

    # 8249
    check_8249 = import_vcfdistance(input = reference, reference = reference,
                sample_name= '8249'
        ).pairwise_distance()

    assert check_8249[3] == 0
    assert all(check_8249[:2] > 0)
