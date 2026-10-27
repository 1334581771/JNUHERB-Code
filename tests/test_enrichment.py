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


class TestTopicRelevance:
    """Relevance-tiered ranking (Core/Related/Review) for enrichment terms."""

    def _kegg_df(self):
        return pd.DataFrame({
            "Term": ["Osteoclast differentiation",
                     "NF-kappa B signaling pathway",
                     "Prostate cancer",
                     "Hepatitis B"],
            "P-value": [1e-6, 1e-9, 1e-4, 1e-3],
            "Category": ["KEGG"] * 4,
            "Genes": ["A/B/C", "A/B/D", "A/B/C", "A/B/C"],
        })

    def test_topic_terms_ranked_core_first(self):
        from jnuherb.enrichment import _prioritise
        out = _prioritise(self._kegg_df(), "osteoarthritis", "kegg")
        # "Osteoclast differentiation" matches topic token "osteo" -> Core,
        # even though NF-kappa B has a smaller P-value.
        assert out.iloc[0]["Term"] == "Osteoclast differentiation"
        assert out.iloc[0]["Relevance"] == "Core"
        assert out.iloc[1]["Relevance"] == "Related"
        # Human-Diseases terms without topic match are pushed to Review.
        review = out[out["Relevance"] == "Review"]
        assert set(review["Term"]) == {"Prostate cancer", "Hepatitis B"}

    def test_brite_classes_assigned(self):
        from jnuherb.enrichment import _prioritise
        out = _prioritise(self._kegg_df(), "osteoarthritis", "kegg")
        classes = dict(zip(out["Term"], out["KEGG Class"]))
        assert classes["Prostate cancer"] == "Human Diseases"
        assert classes["Osteoclast differentiation"] == "Organismal Systems"

    def test_no_topic_keeps_significance_order(self):
        from jnuherb.enrichment import _prioritise
        out = _prioritise(self._kegg_df(), None, "kegg")
        assert set(out["Relevance"]) == {"Related", "Review"}
        # Most significant (NF-kappa B, 1e-9) first among Related.
        related = out[out["Relevance"] == "Related"]
        assert related.iloc[0]["Term"] == "NF-kappa B signaling pathway"

    def test_go_tier_by_topic(self):
        from jnuherb.enrichment import _prioritise
        df = pd.DataFrame({
            "Term": ["inflammatory response", "RNA splicing"],
            "P-value": [1e-3, 1e-9],
            "Category": ["BP", "BP"],
            "Genes": ["A/B", "C/D"],
        })
        out = _prioritise(df, "inflammation", "go")
        assert out.iloc[0]["Term"] == "inflammatory response"
        assert out.iloc[0]["Relevance"] == "Core"


class TestTopTieredFill:
    """Display top-N: Core first, Related fills the remainder, Review last."""

    def test_related_fills_when_core_insufficient(self):
        from jnuherb.enrichment import _prioritise, _top_tiered
        df = pd.DataFrame({
            "Term": ["Osteoclast differentiation", "NF-kappa B signaling pathway",
                     "TNF signaling pathway", "Prostate cancer"],
            "P-value": [1e-4, 1e-9, 1e-8, 1e-3],
            "Category": ["KEGG"] * 4,
            "Genes": ["A/B", "A/C", "A/D", "A/E"],
        })
        tiered = _prioritise(df, "osteoarthritis", "kegg")
        top = _top_tiered(tiered, 10)
        # All 4 terms shown (fewer than 10 available); Core ranks first
        # despite its weaker P-value, Human-Diseases term last.
        assert top["Term"].tolist() == [
            "Osteoclast differentiation",
            "NF-kappa B signaling pathway",
            "TNF signaling pathway",
            "Prostate cancer",
        ]
        assert top.iloc[-1]["Relevance"] == "Review"

    def test_review_only_as_last_resort(self):
        from jnuherb.enrichment import _prioritise, _top_tiered
        related = ["NF-kappa B signaling pathway", "TNF signaling pathway",
                   "IL-17 signaling pathway", "MAPK signaling pathway",
                   "Toll-like receptor signaling pathway"]
        review = ["Prostate cancer", "Breast cancer", "Hepatitis B",
                  "Melanoma", "Colorectal cancer", "Pancreatic cancer"]
        df = pd.DataFrame({
            "Term": ["Osteoclast differentiation"] + related + review,
            # Review terms are deliberately the most significant overall —
            # tiers must still dominate the ordering.
            "P-value": [1e-4] + [1e-8] * 5 + [1e-20] * 6,
            "Category": ["KEGG"] * (6 + 6),
            "Genes": ["A/B"] * (6 + 6),
        })
        tiered = _prioritise(df, "osteoarthritis", "kegg")
        top = _top_tiered(tiered, 10)
        assert len(top) == 10
        assert top.iloc[0]["Term"] == "Osteoclast differentiation"
        # 1 Core + 5 Related fill first; remaining 4 slots go to Review.
        assert (top["Relevance"] == "Review").sum() == 4
        # No Review term outranks a Related term.
        ranks = top["Relevance"].tolist()
        assert ranks.index("Review") > ranks.index("Related")


class TestGeneOverlapDedup:
    """Jaccard-based GO de-replication (enrichplot pairwise_termsim JC)."""

    def test_redundant_terms_collapsed(self):
        from jnuherb.enrichment import _dedup_gene_overlap
        df = pd.DataFrame({
            "Term": ["leukocyte migration", "leukocyte chemotaxis",
                     "RNA splicing"],
            "P-value": [1e-6, 1e-5, 1e-4],
            "Category": ["BP", "BP", "BP"],
            "Genes": ["A/B/C/D", "A/B/C/E", "X/Y/Z"],
        })
        out = _dedup_gene_overlap(df, cutoff=0.7)
        # Jaccard(migration, chemotaxis) = 3/5 = 0.6 < 0.7 -> both kept.
        assert len(out) == 3

    def test_high_overlap_dropped(self):
        from jnuherb.enrichment import _dedup_gene_overlap
        df = pd.DataFrame({
            "Term": ["leukocyte migration", "leukocyte chemotaxis"],
            "P-value": [1e-6, 1e-5],
            "Category": ["BP", "BP"],
            "Genes": ["A/B/C/D/E/F", "A/B/C/D/E/G"],
        })
        out = _dedup_gene_overlap(df, cutoff=0.7)
        # Jaccard = 5/7 ~ 0.71 -> the less significant term is dropped.
        assert out["Term"].tolist() == ["leukocyte migration"]

    def test_categories_independent(self):
        from jnuherb.enrichment import _dedup_gene_overlap
        df = pd.DataFrame({
            "Term": ["t1", "t2"],
            "P-value": [1e-6, 1e-5],
            "Category": ["BP", "CC"],
            "Genes": ["A/B/C/D/E/F", "A/B/C/D/E/G"],
        })
        out = _dedup_gene_overlap(df, cutoff=0.7)
        assert len(out) == 2  # identical gene sets but different ontologies


class TestRunEnrichmentTopic:
    """run_enrichment end-to-end with a topic."""

    def test_kegg_relevance_with_topic(self):
        from jnuherb.enrichment import run_enrichment
        out = run_enrichment(["TNF", "IL1B", "IL6", "PTGS2", "CCL2"],
                             which="kegg", topic="inflammation")
        full = out["kegg"]["full"]
        if not full.empty:
            assert "Relevance" in full.columns
            assert "KEGG Class" in full.columns
            assert full.iloc[0]["Relevance"] == "Core"
