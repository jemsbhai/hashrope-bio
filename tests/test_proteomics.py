"""Tests for proteomics module — synthetic data."""

import pytest
from hashrope import PolynomialHash, Leaf, rope_concat

from hashrope_bio.proteomics.trajectory import frames_identical


@pytest.fixture
def h():
    return PolynomialHash()


class TestFrameComparison:
    def test_identical_frames(self, h):
        frame = b"\x00" * 120  # 10 atoms × 3 floats × 4 bytes
        # Two copies of the same frame
        rope = rope_concat(Leaf(frame, h), Leaf(frame, h), h)
        assert frames_identical(rope, 0, 1, 120, h)

    def test_different_frames(self, h):
        frame_a = b"\x00" * 120
        frame_b = b"\x01" * 120
        rope = rope_concat(Leaf(frame_a, h), Leaf(frame_b, h), h)
        assert not frames_identical(rope, 0, 1, 120, h)
