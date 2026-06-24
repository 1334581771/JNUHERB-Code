#!/usr/bin/env python
"""Example: Query a monomer's targets and trace its TCM/formula origins.

This example demonstrates the step-by-step query API to trace
deoxyshikonin (a naphthoquinone compound) from chemical → protein
targets → containing TCMs → containing formulas.
"""

from __future__ import annotations

import jnuherb

# Step 1: Find the chemical by name
print("=" * 60)
print("Step 1: Find deoxyshikonin in the chemical database")
print("=" * 60)
chem = jnuherb.get_chemicals("Name", "deoxyshikonin")
print(chem[["HVCID", "Name", "CAS_id", "PubChem_id"]].to_string(index=False))

# Step 2: Query its protein targets
print("\n" + "=" * 60)
print("Step 2: Query protein targets of deoxyshikonin")
print("=" * 60)
targets = jnuherb.get_chem_protein_links("HVCID", list(chem["HVCID"]), score=0)
print(f"Found {len(targets)} target interactions")
print(targets.head(10).to_string(index=False))

# Step 3: Get target gene names and descriptions
print("\n" + "=" * 60)
print("Step 3: Get target protein details")
print("=" * 60)
prots = jnuherb.get_proteins("Ensembl_ID", list(targets["Ensembl_ID"]))
print(prots[["Ensembl_ID", "gene_name", "protein_name"]].head(10).to_string(index=False))

# Step 4: Find which TCMs contain this monomer
print("\n" + "=" * 60)
print("Step 4: Find TCMs containing deoxyshikonin")
print("=" * 60)
tcm_links = jnuherb.get_tcm_chem_links("HVCID", list(chem["HVCID"]))
tcms = jnuherb.get_tcm("HVMID", list(tcm_links["HVMID"]))
print(tcms[["HVMID", "cn_name", "pinyin_name"]].to_string(index=False))

# Step 5: Find formulas containing these TCMs
print("\n" + "=" * 60)
print("Step 5: Find formulas containing these TCMs")
print("=" * 60)
formula_links = jnuherb.get_formula_tcm_links("HVMID", list(tcms["HVMID"]))
formulas = jnuherb.get_formula("HVPID", list(formula_links["HVPID"]))
print(f"Found {len(formulas)} formulas")
print(formulas[["HVPID", "name", "Source Document"]].head(5).to_string(index=False))

print("\n✅ Done!")
