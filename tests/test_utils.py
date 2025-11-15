"""Tests for utility functions."""

import pytest
import numpy as np
from utils.common_utils import parse_encode_dict, encode_sequence


class TestParseEncodeDict:
    """Test parse_encode_dict function."""

    def test_default_encoding(self):
        """Test default encoding specification."""
        result = parse_encode_dict(None)
        expected = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        assert result == expected

    def test_list_encoding(self):
        """Test encoding from list."""
        result = parse_encode_dict(['A', 'C', 'G', 'T'])
        expected = {"A": 0, "C": 1, "G": 2, "T": 3}
        assert result == expected

    def test_dict_encoding(self):
        """Test encoding from dict."""
        input_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
        result = parse_encode_dict(input_dict)
        assert result == input_dict

    def test_invalid_type(self):
        """Test invalid type raises TypeError."""
        with pytest.raises(TypeError):
            parse_encode_dict(123)


class TestEncodeSequence:
    """Test encode_sequence function."""

    def test_string_encoding(self):
        """Test encoding from string."""
        sequence = "ACGT"
        result = encode_sequence(sequence)
        assert result.shape == (4, 5)  # 4 bases, 5 categories (A, C, G, N, T)
        assert result.sum() == 4  # One hot encoding

    def test_uppercase_conversion(self):
        """Test that lowercase is converted to uppercase."""
        sequence = "acgt"
        result = encode_sequence(sequence)
        assert result.shape == (4, 5)
        assert result.sum() == 4

    def test_numpy_array_input(self):
        """Test encoding from numpy array."""
        sequence = np.array([b'A', b'C', b'G', b'T'], dtype='|S1')
        result = encode_sequence(sequence)
        assert result.shape == (4, 5)
        assert result.sum() == 4

    def test_ambiguous_bases(self):
        """Test that ambiguous bases are converted to N."""
        sequence = "ACGTN"
        result = encode_sequence(sequence)
        assert result.shape == (5, 5)
        # Check that N is encoded
        assert result[4, 2] == 1  # N should be encoded

    def test_invalid_input_type(self):
        """Test that invalid input raises TypeError."""
        with pytest.raises(TypeError):
            encode_sequence([1, 2, 3, 4])
