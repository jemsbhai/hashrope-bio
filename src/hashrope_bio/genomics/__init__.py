"""Genomics module: FASTA loading, region queries, mutation localization, tandem repeats, gene diff."""

from .fasta import load_fasta_to_rope, load_fasta_chunked, load_fasta_bytes
from .region_query import region_hash, regions_identical
from .mutation import localize_mutation, localize_all_mutations
from .repeats import detect_tandem_repeat, build_repeat_node
from .gene_diff import (
    GeneRegion, Exon, GeneDiff, ExonDiff, GeneDiffReport,
    load_gene_regions, diff_genes, screen_genes,
)
