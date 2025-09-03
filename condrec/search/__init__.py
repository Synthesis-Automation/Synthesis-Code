"""Search utilities for finding similar reactions in a local dataset.

This package currently provides a light-weight in-memory similarity index
based on RDKit Morgan fingerprints. It can be replaced with FAISS or others.
"""

from .similarity import SimilarityIndex, build_index, search_similar  # noqa: F401

