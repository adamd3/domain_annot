"""Unit tests for coordinate_mapper module."""

import pytest
from domain_annot.coordinate_mapper import map_aa_to_genome_coords


def test_positive_strand_mapping():
    # Test positive strand: gene_start=1000, gene_end=1300
    # AA 1..10 -> nt 1000 .. 1029
    g_start, g_end = map_aa_to_genome_coords(1000, 1300, "+", 1, 10)
    assert g_start == 1000
    assert g_end == 1029

    # AA 15..87 for gene starting at 4072028 (+ strand)
    # 4072028 + (15-1)*3 = 4072070
    # 4072028 + (87*3) - 1 = 4072288
    g_start, g_end = map_aa_to_genome_coords(4072028, 4072684, "+", 15, 87)
    assert g_start == 4072070
    assert g_end == 4072288


def test_negative_strand_mapping():
    # Test negative strand: gene_start=1000, gene_end=1300
    # AA 1..10 ->
    # g_start = 1300 - 30 + 1 = 1271
    # g_end = 1300 - 0 = 1300
    g_start, g_end = map_aa_to_genome_coords(1000, 1300, "-", 1, 10)
    assert g_start == 1271
    assert g_end == 1300
    assert g_start <= g_end


def test_invalid_strand():
    g_start, g_end = map_aa_to_genome_coords(1000, 1300, "?", 1, 10)
    assert g_start is None
    assert g_end is None
