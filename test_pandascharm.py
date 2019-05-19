#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

import numpy
import pandas
import dendropy

import Bio.Alphabet
from Bio.AlignIO import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pandas.util.testing import (
    assert_categorical_equal,
    assert_dict_equal,
    assert_frame_equal,
    assert_index_equal,
    assert_produces_warning,
    assert_series_equal)

from pandascharm import (
    frame_as_categorical,
    frame_as_object,
    from_charmatrix,
    to_charmatrix,
    from_bioalignment,
    to_bioalignment,
    from_sequence_dict,
    to_sequence_dict)


class TestAsCategorical():

    frame = pandas.DataFrame({
        't1': ['T', 'G', 'C', 'A', '?'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', 'C', 'A', 'A'],
        't4': ['T', 'G', 'C', 'A', 'A']}, dtype='category')

    def test_unaltered_categories(self):
        assert (
            set(frame_as_categorical(self.frame)['t1'].cat.categories) ==
            set(self.frame['t1'].cat.categories))

    def test_altered_categories(self):
        assert (
            set(frame_as_categorical(self.frame)['t2'].cat.categories) !=
            set(self.frame['t2'].cat.categories))

    def test_add_category(self):
        assert(
            set(
                frame_as_categorical(self.frame, ['-'])['t1'].cat.categories
            ) == {'T', 'G', 'C', 'A', '?', '-'})


class TestAsObject():

    frame_cat = pandas.DataFrame({
        't1': ['T', 'G', 'C', 'A', '?'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', 'C', 'A', 'A'],
        't4': ['T', 'G', 'C', 'A', 'A']}, dtype='category')

    frame_obj = pandas.DataFrame({
        't1': ['T', 'G', 'C', 'A', '?'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', 'C', 'A', 'A'],
        't4': ['T', 'G', 'C', 'A', 'A']}, dtype='object')

    def test_conversion(self):
        assert_frame_equal(frame_as_object(self.frame_cat), self.frame_obj)


class TestCharmatrixConversion():

    dna_charmatrix_string = '3 5\nt1  TCCAA\nt2  TGCAA\nt3  TG-AA\n'
    dna_charmatrix = dendropy.DnaCharacterMatrix.get(
        data=dna_charmatrix_string, schema='phylip')
    dna_frame = pandas.DataFrame({
        't1': ['T', 'C', 'C', 'A', 'A'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', '-', 'A', 'A']}, dtype='category')

    rna_charmatrix_string = '3 5\nt1  UCCAA\nt2  UGCAA\nt3  UG-AA\n'
    rna_charmatrix = dendropy.RnaCharacterMatrix.get(
        data=rna_charmatrix_string, schema='phylip')
    rna_frame = pandas.DataFrame({
        't1': ['U', 'C', 'C', 'A', 'A'],
        't2': ['U', 'G', 'C', 'A', 'A'],
        't3': ['U', 'G', '-', 'A', 'A']}, dtype='category')

    protein_charmatrix_string = '3 5\nt1  VKYPN\nt2  VLYPN\nt3  VL-PN\n'
    protein_charmatrix = dendropy.ProteinCharacterMatrix.get(
        data=protein_charmatrix_string, schema='phylip')
    protein_frame = pandas.DataFrame({
        't1': ['V', 'K', 'Y', 'P', 'N'],
        't2': ['V', 'L', 'Y', 'P', 'N'],
        't3': ['V', 'L', '-', 'P', 'N']}, dtype='category')

    standard_charmatrix_string = '3 5\nt1  01010\nt2  02010\nt3  02-10\n'
    standard_charmatrix = dendropy.StandardCharacterMatrix.get(
        data=standard_charmatrix_string, schema='phylip')
    standard_frame = pandas.DataFrame({
        't1': ['0', '1', '0', '1', '0'],
        't2': ['0', '2', '0', '1', '0'],
        't3': ['0', '2', '-', '1', '0']}, dtype='category')

    def test_from_charmatrix_dna(self):
        assert_frame_equal(
            from_charmatrix(self.dna_charmatrix), self.dna_frame,
            check_categorical=False)

    def test_from_charmatrix_dna_object(self):
        assert_frame_equal(
            from_charmatrix(self.dna_charmatrix, categorical=False),
            frame_as_object(self.dna_frame))

    def test_to_charmatrix_dna(self):
        assert (
            to_charmatrix(self.dna_frame, data_type='dna')
            .as_string('phylip') == self.dna_charmatrix.as_string('phylip'))

    def test_from_charmatrix_rna(self):
        assert_frame_equal(
            from_charmatrix(self.rna_charmatrix), self.rna_frame,
            check_categorical=False)

    def test_to_charmatrix_rna(self):
        assert (
            to_charmatrix(self.rna_frame, data_type='rna')
            .as_string('phylip') == self.rna_charmatrix.as_string('phylip'))

    def test_from_charmatrix_protein(self):
        assert_frame_equal(
            from_charmatrix(self.protein_charmatrix), self.protein_frame,
            check_categorical=False)

    def test_to_charmatrix_protein(self):
        assert (
            to_charmatrix(self.protein_frame, data_type='protein')
            .as_string('phylip') == self.protein_charmatrix
            .as_string('phylip'))

    def test_from_charmatrix_standard(self):
        assert_frame_equal(
            from_charmatrix(self.standard_charmatrix), self.standard_frame,
            check_categorical=False)

    def test_to_charmatrix_standard(self):
        assert (
            to_charmatrix(self.standard_frame, data_type='standard')
            .as_string('phylip') == self.standard_charmatrix
            .as_string('phylip'))

    def test_invalid_data_type(self):
        with pytest.raises(ValueError):
            to_charmatrix(self.standard_frame, data_type='unknown')


class TestBioalignmentConversion():

    def dict_to_bioalignment(d, alphabet='generic_alphabet', sorted=True):
        """
        Create a BioPython MultipleSequenceAlignment
        from a dict with pairs consisting of id and sequence.
        """
        alignment = MultipleSeqAlignment([])
        bio_alphabet = getattr(Bio.Alphabet, alphabet)
        for id, seq in d.items():
            seq_record = SeqRecord(Seq(seq, alphabet=bio_alphabet), id=id)
            alignment.append(seq_record)
        if sorted:
            alignment.sort()
        return alignment

    dna_alignment_dict = {'t1': 'TCCAA', 't2': 'TGCAA', 't3': 'TG-AA'}
    dna_bioalignment = dict_to_bioalignment(
        dna_alignment_dict, alphabet='generic_dna')
    dna_frame = pandas.DataFrame({
        't1': ['T', 'C', 'C', 'A', 'A'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', '-', 'A', 'A']}, dtype='category')

    def test_from_bioalignment_dna(self):
        assert_frame_equal(
            from_bioalignment(self.dna_bioalignment), self.dna_frame)

    def test_to_bioalignment_dna(self):
        assert (
            to_bioalignment(self.dna_frame, alphabet='generic_dna')
            .format('phylip') == self.dna_bioalignment.format('phylip'))

    def test_invalid_alphabet(self):
        with pytest.raises(ValueError):
            to_bioalignment(self.dna_frame, alphabet='dna')


class TestSequenceDictConversion():

    dna_frame = pandas.DataFrame({
        't1': ['T', 'C', 'C', 'A', 'A'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', '-', 'A', 'A']}, dtype='object')

    dna_frame_nan = pandas.DataFrame({
        't1': ['T', 'C', 'C', 'A', 'A'],
        't2': ['T', 'G', 'C', 'A', 'A'],
        't3': ['T', 'G', '-', 'A', numpy.nan]}, dtype='object')

    dna_dict = {'t1': 'TCCAA', 't2': 'TGCAA', 't3': 'TG-AA'}

    def test_from_sequence_dict(self):
        assert_frame_equal(
            from_sequence_dict(self.dna_dict, categorical=False),
            self.dna_frame)

    def test_to_sequence_dict(self):
        assert(to_sequence_dict(self.dna_frame) == self.dna_dict)

    def test_do_sequence_dict_nan(self):
        with pytest.raises(TypeError):
            to_sequence_dict(self.dna_frame_nan)
