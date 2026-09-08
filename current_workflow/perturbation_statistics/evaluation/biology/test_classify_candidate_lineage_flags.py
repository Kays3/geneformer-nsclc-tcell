#!/usr/bin/env python3
"""Regression tests for the identity-only lineage flag classifier."""
from __future__ import annotations

import unittest

import pandas as pd

import classify_candidate_lineage_flags as m


class ClassifyGenesTests(unittest.TestCase):
    def test_named_biology_readme_genes_flag_ambient(self) -> None:
        """The exact genes the biology README calls out by name must flag."""
        named = ["SFTPC", "SFTPB", "NAPSA", "MUC1", "PIGR", "FBLN1", "DCN", "ACKR1"]
        flags = m.classify_genes(pd.Series(named)).set_index("Gene_name")
        for gene in named:
            self.assertEqual(
                flags.loc[gene, "interpretation_class"],
                "4_ambient_doublet_sensitive",
                msg=f"{gene} should flag ambient/doublet-sensitive by identity",
            )

    def test_tcell_anchor_is_pending_not_confirmed(self) -> None:
        flags = m.classify_genes(pd.Series(["GZMK"])).set_index("Gene_name")
        self.assertTrue(
            flags.loc["GZMK", "interpretation_class"].startswith("pending_identity_consistent")
        )
        self.assertNotEqual(flags.loc["GZMK", "interpretation_class"], "1_t_cell_intrinsic_and_donor_stable")

    def test_unrecognized_gene_is_explicitly_pending_not_guessed(self) -> None:
        flags = m.classify_genes(pd.Series(["NOT_A_REAL_GENE_SYMBOL"])).set_index("Gene_name")
        self.assertEqual(
            flags.loc["NOT_A_REAL_GENE_SYMBOL", "interpretation_class"],
            "pending_requires_per_cell_contamination_scoring",
        )

    def test_no_gene_is_silently_dropped(self) -> None:
        genes = ["SFTPC", "GZMK", "SOME_UNKNOWN_GENE", "SFTPC"]
        flags = m.classify_genes(pd.Series(genes))
        self.assertEqual(set(flags["Gene_name"]), {"SFTPC", "GZMK", "SOME_UNKNOWN_GENE"})


if __name__ == "__main__":
    unittest.main()
