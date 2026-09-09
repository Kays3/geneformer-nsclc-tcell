#!/usr/bin/env python3
"""Regression tests for the Priority 1 ambient/lineage cross-expression diagnostic.

Covers the five areas the implementation spec names: count-layer selection,
feature-name mapping, cross-expression score math, the explicit run guard, and
constant-component (doublet/mito) reporting.

The tests build synthetic AnnData objects rather than reaching for the real 592 MiB
h5ad, so they run anywhere. Two of them encode the traps that motivated the spec:

  * ``SiblingObjectRejectionTests`` builds a decoy with the SAME shape as the real
    object and asserts that a shape check passes on it while the loader rejects it.
    A test that only checked dimensions would pass on the wrong file.
  * ``ConstantComponentTests`` asserts the retirement statement flips when handed an
    object with real doublet calls - proving the status is derived from the object,
    not hard-coded prose.
"""
from __future__ import annotations

import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import ambient_risk_diagnostic as m

REAL_CANDIDATES = (
    Path(__file__).resolve().parents[1]
    / "results" / "analysis" / "tables" / "published_gene_donor_support.csv"
)

AMBIENT_ANCHORS = ["SFTPC", "EPCAM", "LYZ"]
TCELL_ANCHORS = ["CD3D", "CD3E", "TRAC"]
OTHER_GENES = ["GENE1", "GENE2"]


def make_object(n_cells=6, symbols=None, with_count_layer=True, with_full_obs=True,
                doublet_values=None, mito_zero=True, counts=None,
                count_layer_name=None, extra_layers=None):
    """A miniature stand-in for the real object, with the same structural contract."""
    symbols = list(symbols if symbols is not None else AMBIENT_ANCHORS + TCELL_ANCHORS + OTHER_GENES)
    n_genes = len(symbols)
    rng = np.random.default_rng(0)
    if counts is None:
        counts = rng.integers(0, 20, size=(n_cells, n_genes)).astype(np.float64)
    counts = np.asarray(counts, dtype=np.float64)

    # X is normalised/log, exactly as on the real object - deliberately NOT the counts.
    library = np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    x = np.log1p(counts / library * 1e4)

    obs = pd.DataFrame({
        "cell_id": [f"cell{i}" for i in range(n_cells)],
        "celltype": ["CD8" if i % 2 else "CD4" for i in range(n_cells)],
        "individual": [f"donor{i % 3}" for i in range(n_cells)],
        "disease": ["LUAD" if i % 2 else "LUSC" for i in range(n_cells)],
        "split": ["test" if i % 2 else "train" for i in range(n_cells)],
    })
    if with_full_obs:
        obs["doublet_status"] = (doublet_values if doublet_values is not None
                                 else ["singlet"] * n_cells)
        obs["pct_counts_mito"] = np.zeros(n_cells) if mito_zero else np.linspace(0.1, 5.0, n_cells)
        obs["total_counts_mito"] = np.zeros(n_cells) if mito_zero else np.arange(n_cells, dtype=float)
        obs["total_counts"] = counts.sum(axis=1)
        obs["n_genes_by_counts"] = (counts > 0).sum(axis=1).astype(float)
    obs.index = [f"cell{i}" for i in range(n_cells)]

    var = pd.DataFrame({"feature_name": symbols, "mito": [False] * n_genes},
                       index=[f"ENSG{i:011d}" for i in range(n_genes)])

    layers = {}
    if with_count_layer:
        layers[count_layer_name or m.COUNT_LAYER] = sp.csr_matrix(counts)
    if extra_layers:
        for name, matrix in extra_layers.items():
            layers[name] = sp.csr_matrix(matrix)
    return ad.AnnData(X=sp.csr_matrix(x), obs=obs, var=var, layers=layers)


class _NoneKeyLayers(dict):
    """Reproduces anndata 0.13.x, whose Layers.keys() yields a literal None beside the
    real layer names. Sorting those keys raises TypeError - a crash that appears only
    against the real object on the real host, never against a locally built AnnData."""


class _Stub:
    """Minimal duck-typed object for contract cases AnnData itself forbids."""
    def __init__(self, var, n_vars, shape=None, layers=None, obs=None):
        self.var, self.n_vars = var, n_vars
        self.shape = shape or (1, n_vars)
        self.layers = layers or {}
        self.obs = obs if obs is not None else pd.DataFrame()


# ---------------------------------------------------------------- 1. count layer

class CountLayerSelectionTests(unittest.TestCase):
    def test_selects_the_count_layer_and_not_x(self) -> None:
        counts = np.array([[5.0, 0.0, 3.0, 1.0, 0.0, 2.0, 4.0, 0.0]])
        adata = make_object(n_cells=1, counts=counts)
        selected = m.select_count_matrix(adata)
        np.testing.assert_allclose(selected.toarray(), counts)
        # The trap: X is present, same shape, numeric - and is NOT what was returned.
        self.assertFalse(np.allclose(np.asarray(adata.X.todense()), counts))

    def test_missing_count_layer_raises_and_does_not_fall_back_to_x(self) -> None:
        adata = make_object(with_count_layer=False)
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.select_count_matrix(adata)
        self.assertIn("layers['count'] is missing", str(ctx.exception))
        self.assertIn("Refusing to fall back to X", str(ctx.exception))

    def test_length_scaled_neighbour_is_named_but_never_substituted(self) -> None:
        """`counts_length_scaled` sits beside `count` on the real object. A fallback
        that quietly used it would report success against the wrong input."""
        adata = make_object(with_count_layer=False,
                            extra_layers={m.LENGTH_SCALED_LAYER: np.ones((6, 8))})
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.select_count_matrix(adata)
        self.assertIn(m.LENGTH_SCALED_LAYER, str(ctx.exception))
        self.assertIn("NOT a substitute", str(ctx.exception))

    def test_non_integral_count_layer_is_rejected(self) -> None:
        adata = make_object()
        adata.layers[m.COUNT_LAYER] = sp.csr_matrix(
            np.asarray(adata.layers[m.COUNT_LAYER].todense()) + 0.5)
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.select_count_matrix(adata)
        self.assertIn("not integral", str(ctx.exception))

    def test_negative_and_non_finite_counts_are_rejected(self) -> None:
        for bad, needle in ((-1.0, "negative"), (np.inf, "non-finite")):
            adata = make_object()
            dense = np.asarray(adata.layers[m.COUNT_LAYER].todense())
            dense[0, 0] = bad
            adata.layers[m.COUNT_LAYER] = sp.csr_matrix(dense)
            with self.assertRaises(m.ObjectContractError) as ctx:
                m.select_count_matrix(adata)
            self.assertIn(needle, str(ctx.exception))


class NoneKeyLayerTests(unittest.TestCase):
    """Regression: anndata 0.13.2 on thinkstation2 exposes a None key in
    layers.keys() while h5py shows only {count, counts_length_scaled}. This crashed
    the first real run with TypeError; the local anndata (0.12.10) does not do it."""

    def test_none_key_beside_a_real_count_layer_does_not_crash(self) -> None:
        counts = np.array([[4.0, 1.0, 0.0]])
        stub = _Stub(var=pd.DataFrame({"feature_name": ["A", "B", "C"]}), n_vars=3,
                     shape=(1, 3),
                     layers=_NoneKeyLayers({m.COUNT_LAYER: sp.csr_matrix(counts),
                                            None: sp.csr_matrix(counts)}))
        np.testing.assert_allclose(m.select_count_matrix(stub).toarray(), counts)

    def test_none_key_alone_raises_the_contract_error_not_typeerror(self) -> None:
        stub = _Stub(var=pd.DataFrame({"feature_name": ["A"]}), n_vars=1, shape=(1, 1),
                     layers=_NoneKeyLayers({None: sp.csr_matrix(np.array([[1.0]]))}))
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.select_count_matrix(stub)
        self.assertIn("layers['count'] is missing", str(ctx.exception))


# ------------------------------------------------------- 2. feature-name mapping

class FeatureMapTests(unittest.TestCase):
    def test_maps_symbols_to_columns_from_the_in_object_var(self) -> None:
        adata = make_object()
        feature_map = m.build_feature_map(adata)
        self.assertEqual(feature_map["SFTPC"], 0)
        self.assertEqual(feature_map["CD3D"], 3)
        self.assertEqual(len(feature_map), adata.n_vars)

    def test_symbols_do_not_match_the_ensembl_var_index(self) -> None:
        """Why the in-object map is needed at all: the var index is Ensembl, so symbol
        lookups against it miss silently rather than erroring."""
        adata = make_object()
        self.assertTrue(all(str(v).startswith("ENSG") for v in adata.var_names))
        self.assertNotIn("CD3D", set(adata.var_names))

    def test_duplicate_symbols_are_rejected_rather_than_silently_collapsed(self) -> None:
        adata = make_object(symbols=["CD3D", "CD3D", "EPCAM", "SFTPC", "LYZ", "CD3E",
                                     "TRAC", "GENE1"])
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.build_feature_map(adata)
        self.assertIn("not unique", str(ctx.exception))
        self.assertIn("CD3D", str(ctx.exception))

    def test_empty_symbol_is_rejected(self) -> None:
        adata = make_object(symbols=["", "EPCAM", "LYZ", "CD3D", "CD3E", "TRAC",
                                     "GENE1", "GENE2"])
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.build_feature_map(adata)
        self.assertIn("empty/NaN", str(ctx.exception))

    def test_length_mismatch_is_rejected(self) -> None:
        stub = _Stub(var=pd.DataFrame({"feature_name": ["A", "B"]}), n_vars=3)
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.build_feature_map(stub)
        self.assertIn("!= n_vars", str(ctx.exception))

    def test_missing_feature_name_column_is_rejected(self) -> None:
        stub = _Stub(var=pd.DataFrame({"other": ["A"]}), n_vars=1)
        with self.assertRaises(m.ObjectContractError):
            m.build_feature_map(stub)


# -------------------------------------------------- 3. cross-expression score math

class ModuleScoreTests(unittest.TestCase):
    def test_score_equals_hand_computed_log1p_cpm_mean(self) -> None:
        counts = np.array([[10.0, 0.0, 5.0, 1.0, 2.0, 0.0, 4.0, 3.0],
                           [0.0, 7.0, 0.0, 6.0, 0.0, 9.0, 1.0, 1.0]])
        adata = make_object(n_cells=2, counts=counts)
        feature_map = m.build_feature_map(adata)
        columns = [feature_map[s] for s in AMBIENT_ANCHORS]

        library = counts.sum(axis=1)
        expected = np.log1p(counts[:, columns] / library[:, None] * m.CPM_SCALE).mean(axis=1)
        actual = m.module_score(sp.csr_matrix(counts), columns)
        np.testing.assert_allclose(actual, expected, rtol=1e-12)

    def test_zero_counts_contribute_zero_not_nan(self) -> None:
        counts = np.zeros((2, 8))
        score = m.module_score(sp.csr_matrix(counts), [0, 1, 2])
        np.testing.assert_allclose(score, np.zeros(2))

    def test_empty_panel_raises_rather_than_returning_nan(self) -> None:
        with self.assertRaises(ValueError):
            m.module_score(sp.csr_matrix(np.ones((2, 3))), [])

    def test_ambient_dominant_cell_scores_above_tcell_dominant_cell(self) -> None:
        """Direction check: the measure must rank a lineage-foreign cell higher."""
        counts = np.zeros((2, 8))
        counts[0, [0, 1, 2]] = 500.0   # ambient anchors dominate cell 0
        counts[1, [3, 4, 5]] = 500.0   # T-cell anchors dominate cell 1
        adata = make_object(n_cells=2, counts=counts)
        counts_matrix, feature_map = m.validate_object(adata)
        frame, provenance = m.compute_cell_scores(adata, counts_matrix, feature_map)
        self.assertGreater(frame.ambient_to_tcell_ratio[0], frame.ambient_to_tcell_ratio[1])
        self.assertGreater(frame.ambient_minus_tcell[0], 0)
        self.assertLess(frame.ambient_minus_tcell[1], 0)
        self.assertEqual(provenance["library_size_source"],
                         "row sums of layers['count'] (not obs['total_counts'])")

    def test_ratio_never_divides_by_zero(self) -> None:
        counts = np.zeros((1, 8))
        counts[0, 0] = 100.0  # ambient only; T-cell module score is exactly zero
        adata = make_object(n_cells=1, counts=counts)
        counts_matrix, feature_map = m.validate_object(adata)
        frame, _ = m.compute_cell_scores(adata, counts_matrix, feature_map)
        self.assertTrue(np.isfinite(frame.ambient_to_tcell_ratio[0]))

    def test_depth_columns_are_labelled_as_depth(self) -> None:
        """total_counts/n_genes_by_counts are sequencing depth, and the column names
        must not let them be read as a contamination measure."""
        adata = make_object()
        counts_matrix, feature_map = m.validate_object(adata)
        frame, _ = m.compute_cell_scores(adata, counts_matrix, feature_map)
        self.assertIn("depth_total_counts", frame.columns)
        self.assertIn("depth_n_genes_by_counts", frame.columns)
        self.assertNotIn("total_counts", frame.columns)

    def test_unresolvable_anchor_panel_raises(self) -> None:
        adata = make_object(symbols=[f"NOSUCH{i}" for i in range(8)])
        counts_matrix = m.select_count_matrix(adata)
        feature_map = m.build_feature_map(adata)
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.compute_cell_scores(adata, counts_matrix, feature_map)
        self.assertIn("did not resolve", str(ctx.exception))


# ------------------------------------------------------------- 4. the run guard

class RunGuardTests(unittest.TestCase):
    def test_refuses_without_the_opt_in_variable(self) -> None:
        with self.assertRaises(SystemExit) as ctx:
            m.require_run_opt_in({})
        self.assertIn(m.ALLOW_RUN_ENV, str(ctx.exception))

    def test_refuses_on_any_value_other_than_one(self) -> None:
        for value in ("0", "true", "yes", ""):
            with self.assertRaises(SystemExit):
                m.require_run_opt_in({m.ALLOW_RUN_ENV: value})

    def test_allows_when_explicitly_opted_in(self) -> None:
        m.require_run_opt_in({m.ALLOW_RUN_ENV: "1"})

    def test_default_path_is_the_confirmed_object(self) -> None:
        self.assertEqual(m.CONFIRMED_H5AD.name,
                         "balanced_lusc_max_7000_per_disease_tcells.h5ad")
        self.assertIn("tcell_luad_lusc_normal_10k_from_atlas", str(m.CONFIRMED_H5AD))


# ------------------------------------------------- 5. constant-component reporting

class ConstantComponentTests(unittest.TestCase):
    def test_states_both_retirements_explicitly(self) -> None:
        report = m.constant_component_report(make_object())
        self.assertEqual(report["doublet_component_status"], "constant_by_construction")
        self.assertEqual(report["doublet_status_values"], ["singlet"])
        self.assertEqual(report["mito_component_status"], "constant_zero_by_construction")

    def test_retirement_is_derived_not_hardcoded(self) -> None:
        """Handed an object with real doublet calls and non-zero mito, the status must
        change. A hard-coded retirement string would pass the test above and fail here."""
        adata = make_object(doublet_values=["singlet", "doublet"] * 3, mito_zero=False)
        report = m.constant_component_report(adata)
        self.assertEqual(report["doublet_component_status"],
                         "present_and_variable__retirement_does_not_apply")
        self.assertEqual(report["doublet_status_values"], ["doublet", "singlet"])
        self.assertEqual(report["mito_component_status"],
                         "not_constant_zero__retirement_does_not_apply")

    def test_report_does_not_present_a_zero_doublet_rate(self) -> None:
        """The failure mode being guarded: a clean-looking 'no doublets, 0% mito'."""
        report = m.constant_component_report(make_object())
        for key in report:
            self.assertNotIn("doublet_rate", key)
            self.assertNotIn("pct_doublets", key)
        self.assertIn("RETIRED", report["retirement_statement"])
        self.assertIn("vacuous", report["retirement_statement"])

    def test_names_depth_columns_as_not_contamination(self) -> None:
        report = m.constant_component_report(make_object())
        self.assertEqual(set(report["depth_columns_not_contamination"]),
                         {"total_counts", "n_genes_by_counts"})
        self.assertIn("NOT contamination", report["depth_note"])

    def test_mito_gene_clause_matches_the_measured_flag_count(self) -> None:
        """On the confirmed object var['mito'] flags ZERO genes, so a statement that it
        "still flags the mitochondrial genes" would be contradicted by the count printed
        beside it. The clause is derived from the count, not asserted alongside it."""
        adata = make_object()  # var['mito'] is all False -> zero flagged
        report = m.constant_component_report(adata)
        self.assertEqual(report["n_var_genes_flagged_mito"], 0)
        self.assertIn("flags ZERO genes", report["retirement_statement"])
        self.assertIn("absent from this object's gene set", report["retirement_statement"])

    def test_mito_gene_clause_flips_when_genes_are_flagged(self) -> None:
        adata = make_object()
        adata.var["mito"] = [True, True] + [False] * (adata.n_vars - 2)
        report = m.constant_component_report(adata)
        self.assertEqual(report["n_var_genes_flagged_mito"], 2)
        self.assertIn("still flags 2 genes", report["retirement_statement"])

    def test_records_that_var_still_flags_mito_genes(self) -> None:
        """So the zeros read as upstream filtering rather than as a bug."""
        report = m.constant_component_report(make_object())
        self.assertTrue(report["var_mito_flag_present"])
        self.assertIsNotNone(report["n_var_genes_flagged_mito"])


# --------------------------------------------- the same-shape sibling object trap

class SiblingObjectRejectionTests(unittest.TestCase):
    """The tokenization input has the IDENTICAL 21,000 x 17,764 shape as the real
    object. Assertions must key on layers/obs, never on dimensions."""

    def _sibling(self):
        return make_object(with_count_layer=False, with_full_obs=False)

    def test_a_shape_check_cannot_tell_the_objects_apart(self) -> None:
        real, sibling = make_object(), self._sibling()
        self.assertEqual(real.shape, sibling.shape)          # a dimension gate passes
        self.assertEqual(real.n_obs, sibling.n_obs)
        self.assertEqual(real.n_vars, sibling.n_vars)

    def test_loader_rejects_the_sibling_on_obs_keys(self) -> None:
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.validate_object(self._sibling())
        message = str(ctx.exception)
        self.assertIn("obs is missing required columns", message)
        self.assertIn("doublet_status", message)
        self.assertIn("dimension check cannot tell the two objects apart", message)

    def test_loader_rejects_a_sibling_that_has_obs_but_no_count_layer(self) -> None:
        with self.assertRaises(m.ObjectContractError) as ctx:
            m.validate_object(make_object(with_count_layer=False))
        self.assertIn("layers['count'] is missing", str(ctx.exception))

    def test_loader_accepts_the_real_contract(self) -> None:
        counts, feature_map = m.validate_object(make_object())
        self.assertEqual(counts.shape, (6, 8))
        self.assertEqual(len(feature_map), 8)


# ------------------------------------------------------------ retained gene level

class LoadCandidatesTests(unittest.TestCase):
    def test_loads_real_local_published_support_table(self) -> None:
        if not REAL_CANDIDATES.exists():
            self.skipTest("published_gene_donor_support.csv not present in this checkout")
        per_gene = m.load_candidates(REAL_CANDIDATES)
        self.assertIn("Gene_name", per_gene.columns)
        self.assertIn("n_comparisons", per_gene.columns)
        sftpc = per_gene.loc[per_gene["Gene_name"] == "SFTPC"]
        self.assertEqual(len(sftpc), 1)
        self.assertTrue(bool(sftpc.iloc[0]["donor_supported_in_all"]))

    def test_missing_column_raises_not_silently_drops(self) -> None:
        import tempfile
        bad = pd.DataFrame({"Gene_name": ["X"], "Ensembl_ID": ["ENSG1"]})
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as handle:
            bad.to_csv(handle.name, index=False)
            path = Path(handle.name)
        try:
            with self.assertRaises(SystemExit):
                m.load_candidates(path)
        finally:
            path.unlink()


class GroupFStatisticTests(unittest.TestCase):
    def test_perfectly_separated_groups_give_large_f(self) -> None:
        counts = sp.csr_matrix(np.array([[0.0], [0.0], [10.0], [10.0]]))
        f = m.group_f_statistic(counts, pd.Series(["a", "a", "b", "b"]))
        self.assertGreater(f[0], 100.0)

    def test_identical_groups_give_near_zero_f(self) -> None:
        counts = sp.csr_matrix(np.array([[5.0], [5.0], [5.0], [5.0]]))
        f = m.group_f_statistic(counts, pd.Series(["a", "a", "b", "b"]))
        self.assertLess(f[0], 1e-6)


if __name__ == "__main__":
    unittest.main()
