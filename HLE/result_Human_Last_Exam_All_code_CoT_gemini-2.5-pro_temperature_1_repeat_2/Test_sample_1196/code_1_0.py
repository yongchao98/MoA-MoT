import numpy as np
from scipy.stats import chi2_contingency

def simulate_gwas_ld():
    """
    Simulates a small GWAS to demonstrate how Linkage Disequilibrium (LD)
    can cause misleading associations.

    Scenario:
    - SNP2 is the true CAUSAL variant for a disease.
    - SNP1 and SNP3 are in high LD (tightly linked) with SNP2 but are NOT causal.
    - SNP4 is on a different chromosome (not in LD) and is NOT causal.
    - We will test the association of each SNP with the disease.
    """
    print("--- Simulation Plan ---")
    print("SNP2 is the true causal SNP for a disease.")
    print("SNP1 and SNP3 are non-causal but tightly linked to SNP2.")
    print("SNP4 is non-causal and on a different chromosome (unlinked).")
    print("We expect SNP1, SNP2, and SNP3 to show a significant association,")
    print("making the association for SNP1 and SNP3 'misleading'.\n")

    # --- 1. Define Haplotypes and Frequencies ---
    # Haplotype format: [Allele_SNP1, Allele_SNP2_Causal, Allele_SNP3]
    # 'T' at SNP2 is the risk allele.
    haplotypes = {
        'H1': {'alleles': ['A', 'G', 'C'], 'freq': 0.7}, # Common, non-risk
        'H2': {'alleles': ['G', 'T', 'T'], 'freq': 0.3}  # Less common, risk
    }
    # Alleles for unlinked SNP4
    snp4_alleles = {'A': 0.5, 'G': 0.5}

    # --- 2. Generate a Population ---
    n_individuals = 2000
    population = []
    haplotype_names = list(haplotypes.keys())
    haplotype_freqs = [h['freq'] for h in haplotypes.values()]

    for _ in range(n_individuals):
        # Assign two haplotypes for the linked region (Chr 1)
        h1_idx, h2_idx = np.random.choice(len(haplotype_names), 2, p=haplotype_freqs)
        geno_snp1 = haplotypes[haplotype_names[h1_idx]]['alleles'][0] + haplotypes[haplotype_names[h2_idx]]['alleles'][0]
        geno_snp2_causal = haplotypes[haplotype_names[h1_idx]]['alleles'][1] + haplotypes[haplotype_names[h2_idx]]['alleles'][1]
        geno_snp3 = haplotypes[haplotype_names[h1_idx]]['alleles'][2] + haplotypes[haplotype_names[h2_idx]]['alleles'][2]

        # Assign genotype for unlinked SNP4 (Chr 2)
        a1_snp4, a2_snp4 = np.random.choice(list(snp4_alleles.keys()), 2, p=list(snp4_alleles.values()))
        geno_snp4 = a1_snp4 + a2_snp4

        # --- 3. Assign Phenotype (Disease Status) based on Causal SNP2 ---
        risk_allele_count = geno_snp2_causal.count('T')
        disease_prob = 0.10 # Base risk
        if risk_allele_count == 1:
            disease_prob = 0.25 # Heterozygous risk
        elif risk_allele_count == 2:
            disease_prob = 0.50 # Homozygous risk

        disease_status = np.random.rand() < disease_prob

        population.append({
            'snp1': geno_snp1, 'snp2_causal': geno_snp2_causal,
            'snp3': geno_snp3, 'snp4': geno_snp4,
            'diseased': disease_status
        })

    # --- 4. Perform Association Test for each SNP ---
    print("--- Association Test Results (Chi-Squared p-values) ---")
    for snp_name in ['snp1', 'snp2_causal', 'snp3', 'snp4']:
        # We test a dominant model: presence of risk-associated allele vs disease
        # For SNP1 and SNP3, the risk-associated allele is 'G' and 'T' respectively,
        # as they are on the same haplotype as the causal 'T' at SNP2.
        # For SNP4, neither is associated, we'll just test 'G'.
        risk_allele_map = {'snp1': 'G', 'snp2_causal': 'T', 'snp3': 'T', 'snp4': 'G'}
        risk_allele = risk_allele_map[snp_name]

        # Create 2x2 contingency table:
        #           [ Diseased | Healthy ]
        # Has Risk   [          |         ]
        # No Risk    [          |         ]
        has_risk_diseased = len([p for p in population if risk_allele in p[snp_name] and p['diseased']])
        has_risk_healthy = len([p for p in population if risk_allele in p[snp_name] and not p['diseased']])
        no_risk_diseased = len([p for p in population if risk_allele not in p[snp_name] and p['diseased']])
        no_risk_healthy = len([p for p in population if risk_allele not in p[snp_name] and not p['diseased']])

        contingency_table = np.array([
            [has_risk_diseased, has_risk_healthy],
            [no_risk_diseased, no_risk_healthy]
        ])

        chi2, p_val, _, _ = chi2_contingency(contingency_table)
        print(f"Association for {snp_name.upper()}: p-value = {p_val:.4e}")

    print("\n--- Conclusion ---")
    print("As shown, the non-causal SNPs (SNP1, SNP3) that are tightly linked to the")
    print("causal SNP (SNP2) show a highly significant p-value, similar to the causal SNP.")
    print("This is a 'misleading association'. The unlinked SNP4 shows no association.")
    print("This demonstrates why option A is the most likely scenario to cause this issue.")


if __name__ == '__main__':
    simulate_gwas_ld()