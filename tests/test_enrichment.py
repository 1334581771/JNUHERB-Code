"""Tests for enrichment result post-processing."""

from __future__ import annotations

import pandas as pd

from jnuherb.enrichment import _tidy_enr_result


def _enr_df(terms: list[str]) -> pd.DataFrame:
    return pd.DataFrame({
        "Term": terms,
        "Overlap": ["3/10"] * len(terms),
        "P-value": [1e-5] * len(terms),
        "Adjusted P-value": [1e-4] * len(terms),
    })


class TestDiseaseTermDisplay:
    def test_id_suffix_stripped(self):
        df = _enr_df(["acute kidney injury (MONDO:0002492)",
                      "Fatty Liver (MESH:D005234)",
                      "multiple sclerosis, susceptibility to (OMIM:612594)"])
        out = _tidy_enr_result(df, ["A", "B"], "Disease")
        assert out["Term"].tolist() == [
            "Acute kidney injury",
            "Fatty Liver",
            "Multiple sclerosis, susceptibility to",
        ]

    def test_first_letter_capitalized_only(self):
        """Sentence case: first letter up, the rest of the string untouched."""
        df = _enr_df(["non-alcoholic Fatty Liver Disease (MESH:D065626)",
                      "proteostasis deficiencies (MONDO:0021179)"])
        out = _tidy_enr_result(df, ["A"], "Disease")
        assert out["Term"].tolist() == [
            "Non-alcoholic Fatty Liver Disease",
            "Proteostasis deficiencies",
        ]

    def test_other_categories_untouched(self):
        """KEGG/GO term conventions must not be rewritten."""
        df = _enr_df(["ABC transporters", "xenobiotic transport (GO:0042908)"])
        out = _tidy_enr_result(df, ["A"], "KEGG")
        assert out["Term"].tolist() == ["ABC transporters",
                                        "xenobiotic transport (GO:0042908)"]

    def test_empty_result_survives(self):
        assert _tidy_enr_result([], ["A"], "Disease").empty
