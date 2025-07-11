# A pedagogical script to illustrate why tightly linked SNPs can be misleading in GWAS

# Let's define our SNPs and their properties
# We will simulate a scenario with 5 SNPs.
# SNP_C is the true CAUSAL SNP for a trait.
# SNP_L1, SNP_L2, and SNP_L3 are tightly LINKED to SNP_C within the same LD block.
# SNP_U is UNLINKED (e.g., on another chromosome).

snps = {
    "SNP_L1": {"chromosome": 1, "position": 10010, "is_causal": False},
    "SNP_C":  {"chromosome": 1, "position": 10050, "is_causal": True},
    "SNP_L2": {"chromosome": 1, "position": 10085, "is_causal": False},
    "SNP_L3": {"chromosome": 1, "position": 10120, "is_causal": False},
    "SNP_U":  {"chromosome": 2, "position": 50400, "is_causal": False}
}

# Linkage Disequilibrium (LD) is often measured by r^2.
# r^2 = 1 means perfect correlation (the SNPs are always inherited together).
# r^2 = 0 means no correlation.
# Let's define the LD (as r^2) of each SNP with the causal SNP_C.
ld_with_causal = {
    "SNP_L1": 0.95, # Tightly linked
    "SNP_C":  1.0,  # A SNP is perfectly correlated with itself
    "SNP_L2": 0.98, # Tightly linked
    "SNP_L3": 0.90, # Tightly linked
    "SNP_U":  0.01  # Essentially unlinked
}

# Now, let's simulate the results from a Genome-Wide Association Study (GWAS).
# In a GWAS, we calculate a p-value for the association of each SNP with the trait.
# A very small p-value (e.g., < 5e-8) suggests a significant association.
# The true causal SNP (SNP_C) will have a highly significant p-value.
p_value_causal = 1e-10

# The p-values for the linked SNPs will also be highly significant because they are
# correlated with (i.e., in high LD with) the causal SNP. They "hitchhike" on its signal.
# The unlinked SNP will not have a significant p-value.
gwas_results = {}
for snp_name in snps:
    # A simplified model: p-value gets less significant as LD decreases.
    # This is just for illustration. Real p-value relationships are more complex.
    ld_r2 = ld_with_causal[snp_name]
    if ld_r2 > 0.1:
      # If linked, the p-value is a function of the causal p-value and the LD.
      # A simple way to model this is p_observed = p_causal / r^2
      p_observed = p_value_causal / ld_r2
    else:
      # If unlinked, p-value is random and non-significant.
      p_observed = 0.45

    gwas_results[snp_name] = p_observed

# --- Output the Explanation and Results ---
print("--- GWAS Simulation: The Misleading Nature of Linkage Disequilibrium ---")
print("\nContext: A complex trait is influenced by a single Causal SNP (SNP_C).")
print("However, SNP_C is located in an LD Block, tightly linked with three other non-causal SNPs (SNP_L1, SNP_L2, SNP_L3).")
print("This scenario corresponds to Answer Choice A.\n")

print("GWAS Results (P-values):")
print("--------------------------")
for snp_name, p_val in gwas_results.items():
    causal_status = "(CAUSAL)" if snps[snp_name]['is_causal'] else "(Non-Causal)"
    print(f"{snp_name:<7} (Chr {snps[snp_name]['chromosome']}):\tp-value = {p_val:.2e} {causal_status}")
print("--------------------------\n")

print("Analysis:")
print("Observe that not only the true CAUSAL SNP (SNP_C) shows a highly significant p-value,")
print("but all the tightly linked non-causal SNPs (SNP_L1, SNP_L2, SNP_L3) ALSO show highly significant p-values.")
print("\nThis creates a misleading situation because:")
print("1. It's impossible to tell from this data alone which of the four significant SNPs is the true cause.")
print("   The association signal is 'smeared' across the entire LD block.")
print("2. It creates a very strong, concentrated signal at a single genetic locus. If the true genetic basis")
print("   of the trait involves many weaker signals from other genomic regions, this single strong hit might")
print("   mislead researchers into overestimating its importance and ignoring other true, but weaker, associations.")

print("\nConclusion:")
print("Therefore, a combination of multiple, tightly linked SNPs within a single LD block (Choice A)")
print("is the scenario most likely to provide a misleading association due to linkage disequilibrium.")