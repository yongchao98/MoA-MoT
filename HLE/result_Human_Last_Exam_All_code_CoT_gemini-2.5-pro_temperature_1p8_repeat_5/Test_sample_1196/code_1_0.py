import numpy as np
from scipy import stats

def simulate_gwas_scenario():
    """
    Simulates a GWAS scenario to demonstrate the effect of Linkage Disequilibrium (LD).

    This simulation shows why 'Three SNPs tightly linked within a single LD block'
    is most likely to provide a misleading association.
    """
    n_individuals = 2000
    # Minor allele frequency for the causal SNP
    causal_maf = 0.25

    # --- Genotype Simulation ---
    # 1. Create the true causal SNP
    # Genotypes are coded as 0, 1, or 2 copies of the minor allele.
    causal_snp = np.random.binomial(2, causal_maf, n_individuals)

    # 2. Create three SNPs in high LD with the causal SNP (Scenario A)
    # We simulate high LD by copying the causal SNP's genotypes and adding minimal 'noise'.
    # A low flip_rate means high LD.
    flip_rate = 0.01  # 1% chance to not match the causal allele, creating >98% LD
    ld_snp_1 = np.abs(causal_snp - np.random.binomial(2, flip_rate, n_individuals)) % 3
    ld_snp_2 = np.abs(causal_snp - np.random.binomial(2, flip_rate, n_individuals)) % 3
    ld_snp_3 = np.abs(causal_snp - np.random.binomial(2, flip_rate, n_individuals)) % 3

    # 3. Create a SNP not in LD (e.g., on another chromosome, Scenario C)
    # This SNP is generated independently.
    unlinked_maf = 0.3
    unlinked_snp = np.random.binomial(2, unlinked_maf, n_individuals)

    # --- Trait Simulation ---
    # The trait is a function of the causal SNP's genotype plus random noise.
    # Effect size of the causal SNP
    effect_size = 0.4
    noise = np.random.normal(0, 1, n_individuals)
    trait = effect_size * causal_snp + noise

    # --- Association Testing ---
    snps_to_test = {
        "Causal SNP (True Signal)": causal_snp,
        "Linked SNP 1 (Misleading Signal)": ld_snp_1,
        "Linked SNP 2 (Misleading Signal)": ld_snp_2,
        "Linked SNP 3 (Misleading Signal)": ld_snp_3,
        "Unlinked SNP (No Signal)": unlinked_snp
    }

    print("--- GWAS Association Test Simulation ---")
    print(f"Simulating {n_individuals} individuals.")
    print("A lower p-value indicates a stronger association with the trait.\n")

    for name, genotypes in snps_to_test.items():
        # Perform linear regression: trait ~ genotype
        # This is a simplified GWAS association test.
        lin_reg = stats.linregress(x=genotypes, y=trait)
        
        # The equation is trait = intercept + slope * genotype
        slope = lin_reg.slope
        intercept = lin_reg.intercept
        p_value = lin_reg.pvalue

        print(f"Testing SNP: {name}")
        print(f"    Association Equation: trait = {intercept:.4f} + ({slope:.4f} * genotype)")
        print(f"    P-value: {p_value:.2e}")
        # Add interpretation based on p-value
        if p_value < 1e-5:
            print("    Result: Strong statistical association found.\n")
        else:
            print("    Result: No significant association found.\n")

    print("--- Conclusion ---")
    print("As shown, the three SNPs in high LD with the Causal SNP also show a")
    print("strong statistical association. If the Causal SNP was not measured,")
    print("the signals from the Linked SNPs would be considered 'misleading' because")
    print("they are only proxies for the true signal, not the cause themselves.")
    print("This demonstrates the problem described in choice A.")


simulate_gwas_scenario()
<<<A>>>