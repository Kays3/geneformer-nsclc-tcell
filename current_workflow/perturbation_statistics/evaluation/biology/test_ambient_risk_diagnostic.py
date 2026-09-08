#!/usr/bin/env python3
"""Tests for the parts of ambient_risk_diagnostic.py that do NOT need the
(currently unreachable) h5ad: candidate loading, symbol-map loading, and the
group F-statistic helper. main()'s h5ad-dependent feature computation is
untested - it cannot be, without the source single-cell object."""
from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import ambient_risk_diagnostic as m

REAL_CANDIDATES = (
    Path(__file__).resolve().parents[1]
    / "results" / "analysis" / "tables" / "published_gene_donor_support.csv"
)
REAL_GENEFORMER_DIR = Path.home() / "Geneformer" / "geneformer"


class LoadCandidatesTests(unittest.TestCase):
    def test_loads_real_local_published_support_table(self) -> None:
        if not REAL_CANDIDATES.exists():
            self.skipTest("published_gene_donor_support.csv not present in this checkout")
        per_gene = m.load_candidates(REAL_CANDIDATES)
        self.assertIn("Gene_name", per_gene.columns)
        self.assertIn("n_comparisons", per_gene.columns)
        # SFTPC is donor-supported in every one of the (up to 6) comparisons it appears in.
        sftpc = per_gene.loc[per_gene["Gene_name"] == "SFTPC"]
        self.assertEqual(len(sftpc), 1)
        self.assertTrue(bool(sftpc.iloc[0]["donor_supported_in_all"]))

    def test_missing_column_raises_not_silently_drops(self) -> None:
        bad = pd.DataFrame({"Gene_name": ["X"], "Ensembl_ID": ["ENSG1"]})
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as f:
            bad.to_csv(f.name, index=False)
            path = Path(f.name)
        try:
            with self.assertRaises(SystemExit):
                m.load_candidates(path)
        finally:
            path.unlink()


class LoadSymbolMapsTests(unittest.TestCase):
    def test_loads_real_local_gc104m_dictionary(self) -> None:
        if not (REAL_GENEFORMER_DIR / "gene_name_id_dict_gc104M.pkl").exists():
            self.skipTest("gc104M dictionaries not present on this machine")
        symbol_to_ensembl, ensembl_to_symbol = m.load_symbol_maps(REAL_GENEFORMER_DIR)
        self.assertGreater(len(symbol_to_ensembl), 1000)
        self.assertIn("SFTPC", symbol_to_ensembl)
        ensembl_id = symbol_to_ensembl["SFTPC"]
        self.assertEqual(ensembl_to_symbol[ensembl_id], "SFTPC")

    def test_missing_dictionary_raises(self) -> None:
        with self.assertRaises(SystemExit):
            m.load_symbol_maps(Path("/nonexistent/path/for/testing"))


class GroupFStatisticTests(unittest.TestCase):
    def test_perfectly_separated_groups_give_large_f(self) -> None:
        counts = sp.csr_matrix(np.array([[0.0], [0.0], [10.0], [10.0]]))
        labels = pd.Series(["a", "a", "b", "b"])
        f = m.group_f_statistic(counts, labels)
        self.assertGreater(f[0], 100.0)

    def test_identical_groups_give_near_zero_f(self) -> None:
        counts = sp.csr_matrix(np.array([[5.0], [5.0], [5.0], [5.0]]))
        labels = pd.Series(["a", "a", "b", "b"])
        f = m.group_f_statistic(counts, labels)
        self.assertLess(f[0], 1e-6)


class MainGuardTests(unittest.TestCase):
    def test_main_refuses_to_run_against_unverified_path(self) -> None:
        """The h5ad-dependent path must not silently execute against the
        unverified default H5AD path - this is the guard the module docstring
        promises, not a placeholder."""
        with self.assertRaises(SystemExit) as ctx:
            m.main()
        self.assertIn("unrun", str(ctx.exception).lower())


if __name__ == "__main__":
    unittest.main()
