import numpy as np
import pytest
from warnings import warn

from methlab.ibdpainting import *

input = 'tests/test_data/panel_to_test.vcf.gz'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.vcf.gz'
chr1 = 'tests/test_data/reference_panel_chr1.vcf.gz'

input = 'tests/test_data/panel_to_test.vcf.gz'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.vcf.gz'
chr1 = 'tests/test_data/reference_panel_chr1.vcf.gz'


out = load_genotypes(
    input = input,
    reference = reference,
    sample_name = 'S2.15.002'
)

windows = out.split_into_windows(10000)

windows['Chr2:80000-90000'].pairwise_distance()



def test_load_genotypes_gives_right_output():
    """
    Check that load_genotypes gives the right answers when it should.
    """
    out = load_genotypes(
        input = input,
        reference = reference,
        sample_name = 'S2.15.002'
    )
    real_names = np.array(['S2.15.002', '1158', '6024', '6184', '8249'])
    assert all(
        [ x == y for x,y in zip(out.samples, real_names) ]
    )

    assert len(out.chr) == 550
    assert len(out.pos) == 550
    assert out.geno.shape == (550, 5, 2)

def test_load_genotypes_fails_if_missing_sample():
    """
    Check that load_genotypes fails if the sample name is not in the input VCF
    """
    with pytest.raises(Exception):
        load_genotypes(
            input = input,
            reference = reference,
            sample_name = 'not_a_real_sample_name'
        )

def test_load_genotypes_fails_if_contigs_dont_match():
    with pytest.raises(Exception):
        load_genotypes(
            input = input,
            reference = chr1,
            sample_name = '1158_2'
        )

def test_split_into_windows_functions():
    vcfd = load_genotypes(
            input = input,
            reference = reference,
            sample_name = 'S2.15.002'
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
    check_1158 = load_genotypes(
        input = ref_vcf,
        reference = reference,
        sample_name= '1158'
        ).pairwise_distance()

    assert check_1158[0] == 0
    assert all(check_1158[1:] > 0)

    # 6024
    check_6024 = load_genotypes(input = ref_vcf, reference = reference,
                sample_name= '6024'
        ).pairwise_distance()

    assert check_6024[1] == 0
    assert all(check_6024[[0,2,3]] > 0)

    # 6184
    check_6184 = load_genotypes(input = ref_vcf, reference = reference,
                sample_name= '6184'
        ).pairwise_distance()

    assert check_6184[2] == 0
    assert all(check_6184[[0,1,3]] > 0)

    # 8249
    check_8249 = load_genotypes(input = ref_vcf, reference = reference,
                sample_name= '8249'
        ).pairwise_distance()

    assert check_8249[3] == 0
    assert all(check_8249[:2] > 0)

def test_ibdpainting():
    ibd = ibd_table(
        input=ref_vcf,
        reference=reference,
        sample_name='1158',
        window_size=1000
    )
    # Check the dataframe is the right shape
    assert ibd.shape == (200, 5)
    # Check that the column for the true parent is all zeroes or -9
    assert all(
        (ibd['1158'] == 0) | (ibd['1158'] == -9)
    )
    # Check that a non-parent are not all -9.
    assert any(ibd['8249'] != 0)

"""
methlab ibdpainting \
    --input tests/test_data/reference_panel.vcf.gz \
    --reference tests/test_data/reference_panel.vcf.gz \
    --sample_name 1158 \
    --window_size 1000 \
    --outdir tests/test_output \
    --expected_match 1158 \
    --keep_ibd_table
"""