"""Genomics module: FASTA loading, region queries, mutation localization, tandem repeats."""

from .fasta import load_fasta_to_rope, load_fasta_chunked
from .region_query import region_hash, regions_identical
from .mutation import localize_mutation, localize_all_mutations
from .repeats import detect_tandem_repeat, build_repeat_node
