#!/usr/bin/env python3
"""Regression tests for the P0 gate-scoring bug: a hardcoded literal standing
in for a column that every other branch computes.

evaluate_donor_gene_effects.py:539 used to read a literal `"n_ci_excludes_zero": 0`
inside the `published_leading_evaluable` branch of score_gates(), while the
identical field was computed from data in every other scored_set via _set_row().
The published branch silently reported zero CI-excluding genes regardless of
what the underlying effects actually showed."""
from __future__ import annotations

import unittest

import pandas as pd

import evaluate_donor_gene_effects as m


def _effects_row(comparison_label, source, goal, gene, n_donors, ci_excludes_zero, passes):
    return {
        "comparison_label": comparison_label,
        "source": source,
        "goal": goal,
        "Gene_name": gene,
        "Ensembl_ID": f"ENSG_{gene}",
        "n_donors_evaluable": n_donors,
        "min_donors_required": m.min_donors_required(source),
        "donor_supported": n_donors >= m.min_donors_required(source),
        "n_cells_total": 100,
        "median_donor_effect": 0.1 if passes else -0.1,
        "sign_concordance": 1.0 if passes else 0.0,
        "ci_excludes_zero": ci_excludes_zero,
        "same_direction_pass_70": passes,
    }


class PublishedSupportCarriesCiExcludesZero(unittest.TestCase):
    def setUp(self) -> None:
        # LUSC has a 3-donor floor; every gene here clears it.
        self.effects = pd.DataFrame(
            [
                _effects_row("LUSC → LUAD", "lusc", "luad", "GENE_A", 4, True, True),
                _effects_row("LUSC → LUAD", "lusc", "luad", "GENE_B", 4, False, True),
                _effects_row("LUSC → LUAD", "lusc", "luad", "GENE_C", 4, True, True),
            ]
        )
        self.top_genes = pd.DataFrame(
            [
                {"comparison_label": "LUSC → LUAD", "Gene_name": "GENE_A", "Ensembl_ID": "ENSG_GENE_A"},
                {"comparison_label": "LUSC → LUAD", "Gene_name": "GENE_B", "Ensembl_ID": "ENSG_GENE_B"},
                {"comparison_label": "LUSC → LUAD", "Gene_name": "GENE_C", "Ensembl_ID": "ENSG_GENE_C"},
            ]
        )

    def test_published_support_reports_actual_ci_excludes_zero(self) -> None:
        published = m.published_support(self.effects, self.top_genes)
        self.assertIn("ci_excludes_zero", published.columns)
        by_gene = published.set_index("Gene_name")["ci_excludes_zero"].to_dict()
        self.assertEqual(by_gene, {"GENE_A": True, "GENE_B": False, "GENE_C": True})

    def test_missing_published_gene_reports_false_not_a_crash(self) -> None:
        top_genes = pd.concat(
            [
                self.top_genes,
                pd.DataFrame([{"comparison_label": "LUSC → LUAD", "Gene_name": "NOT_PRESENT", "Ensembl_ID": ""}]),
            ],
            ignore_index=True,
        )
        published = m.published_support(self.effects, top_genes)
        row = published.loc[published["Gene_name"] == "NOT_PRESENT"].iloc[0]
        self.assertFalse(row["ci_excludes_zero"])
        self.assertFalse(row["evaluable_for_direction_gate"])

    def test_score_gates_published_branch_computes_not_hardcodes(self) -> None:
        """The regression test: two of three evaluable published genes have
        ci_excludes_zero=True. A hardcoded literal would report 0 here."""
        lodo = pd.DataFrame(
            columns=["ranking", "comparison_label", "level", "gene_retained_in_70pct_folds"]
        )
        published = m.published_support(self.effects, self.top_genes)
        gates = m.score_gates(self.effects, lodo, self.top_genes, published)
        row = gates["by_comparison"]
        pub_row = row.loc[row["scored_set"] == "published_leading_evaluable"].iloc[0]
        self.assertEqual(pub_row["n_ci_excludes_zero"], 2)

    def test_no_scored_set_hardcodes_a_field_set_row_computes(self) -> None:
        """Generalisable guard: every row score_gates() emits for a comparison
        must come from _set_row(), so no branch can drift into a literal for a
        field the others compute from data. This walks the actual gate table
        rather than re-checking one field by name."""
        lodo = pd.DataFrame(
            columns=["ranking", "comparison_label", "level", "gene_retained_in_70pct_folds"]
        )
        published = m.published_support(self.effects, self.top_genes)
        gates = m.score_gates(self.effects, lodo, self.top_genes, published)
        by_comp = gates["by_comparison"]
        expected_fields = set(
            m._set_row("x", "y", self.effects.iloc[0:0], "lusc", 3).keys()
        )
        for row in by_comp.to_dict(orient="records"):
            present = expected_fields & set(row.keys())
            self.assertEqual(present, expected_fields)


class FailMarginGenesTests(unittest.TestCase):
    def test_at_threshold_pass_has_one_gene_margin(self) -> None:
        # 35/50 = 0.70 exactly: losing one gene drops it to 34/50 = 0.68, below the gate.
        self.assertEqual(m._fail_margin_genes(50, 35, True), 1.0)

    def test_two_gene_margin_above_threshold(self) -> None:
        # 36/50 = 0.72: losing two genes drops it to 34/50 = 0.68, below the gate.
        self.assertEqual(m._fail_margin_genes(50, 36, True), 2.0)

    def test_failing_gate_has_no_margin(self) -> None:
        import math

        self.assertTrue(math.isnan(m._fail_margin_genes(50, 30, False)))

    def test_empty_set_has_no_margin(self) -> None:
        import math

        self.assertTrue(math.isnan(m._fail_margin_genes(0, 0, False)))


if __name__ == "__main__":
    unittest.main()
