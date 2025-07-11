import math

def demonstrate_pgs_limitations():
    """
    This function demonstrates why the predictive ability of a Polygenic Score (PGS)
    is necessarily lower than the SNP heritability of a trait.
    """
    # --- Parameters ---
    # SNP heritability (h²_snp): The theoretical maximum variance explained by all SNPs.
    # This is a property of the trait itself.
    h2_snp = 0.50

    # GWAS sample size (N): The number of individuals in the study used to estimate SNP effects.
    N = 150000

    # Effective number of causal variants (M): A measure of the genetic complexity of the trait.
    M = 10000

    print("This script illustrates the relationship between SNP heritability and Polygenic Score (PGS) predictive ability.")
    print("-" * 70)
    print(f"1. SNP Heritability (h²_snp): The theoretical maximum predictive power.")
    print(f"   - Assumed value: {h2_snp}")
    print("\n2. Polygenic Score (PGS): A practical predictor built from a GWAS.")
    print(f"   - To build a PGS, we estimate SNP effects from a GWAS of a finite size.")
    print(f"   - Assumed GWAS sample size (N): {N:,}")
    print(f"   - Assumed number of effective causal variants (M): {M:,}")
    print("-" * 70)

    # --- Calculation ---
    # The predictive ability (R²) of a PGS is limited by estimation error from the finite GWAS sample.
    # The expected maximum R² can be approximated by the following formula:
    # R²_pgs = h²_snp / (1 + M / (N * h²_snp))
    # The term `M / (N * h²_snp)` represents the ratio of noise to signal.
    # As long as N is finite, this ratio is > 0, making the denominator > 1.
    
    if N * h2_snp == 0:
        expected_r2_pgs = 0
    else:
        expected_r2_pgs = h2_snp / (1 + M / (N * h2_snp))

    # --- Output ---
    print("The formula for expected maximum PGS predictive ability (R²) is:")
    print("R² = h²_snp / (1 + M / (N * h²_snp))")
    print("\nPlugging in our assumed values, the equation becomes:")
    print(f"R² = {h2_snp} / (1 + {M} / ({N} * {h2_snp}))")
    
    intermediate_calc = M / (N * h2_snp)
    print(f"R² = {h2_snp} / (1 + {intermediate_calc:.4f})")

    denominator = 1 + intermediate_calc
    print(f"R² = {h2_snp} / {denominator:.4f}")
    
    print(f"\nFinal Calculated PGS R² = {expected_r2_pgs:.4f}")

    print("\n--- Conclusion ---")
    print(f"The theoretical maximum (SNP Heritability = {h2_snp}) is higher than the")
    print(f"practical predictive ability (PGS R² = {expected_r2_pgs:.4f}).")
    print("\nThis demonstrates that because a PGS relies on *estimated* effects from a finite sample, its predictive power is necessarily lower than the true SNP heritability.")

if __name__ == '__main__':
    demonstrate_pgs_limitations()
