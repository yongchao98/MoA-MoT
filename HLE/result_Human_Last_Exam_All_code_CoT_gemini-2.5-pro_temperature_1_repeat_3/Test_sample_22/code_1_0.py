def demonstrate_pgs_prediction_limit():
    """
    Demonstrates that the predictive ability of a Polygenic Score (R^2)
    is lower than the total SNP heritability (h2_snp).
    """

    # --- Input Parameters ---
    # SNP heritability: The total proportion of variance explained by all SNPs.
    # Let's assume a trait with 50% heritability.
    h2_snp = 0.5

    # GWAS sample size: The number of individuals in the discovery study.
    # A large, but finite, sample size.
    N = 150000

    # Effective number of independent SNPs. This is typically in the tens of thousands
    # for European ancestry populations.
    M = 60000

    # --- Calculation ---
    # The formula for the expected R^2 of a PGS is:
    # E[R^2] = h2_snp / (1 + M / (N * h2_snp))
    # This formula shows how the maximum possible R^2 (which is h2_snp) is
    # reduced by a factor related to GWAS sample size and number of SNPs.

    numerator = h2_snp
    denominator = 1 + (M / (N * h2_snp))
    expected_r2_pgs = numerator / denominator

    # --- Output ---
    print(f"Given Parameters:")
    print(f"  - Total SNP Heritability (h2_snp): {h2_snp}")
    print(f"  - GWAS Sample Size (N): {N}")
    print(f"  - Effective Number of SNPs (M): {M}")
    print("\n" + "="*40 + "\n")
    print(f"Calculation of Expected PGS Predictive Ability (R^2):")
    print(f"  - Numerator (h2_snp): {numerator}")
    print(f"  - Denominator (1 + M / (N * h2_snp)): {denominator:.4f}")
    print(f"  - Expected PGS R^2 = {numerator} / {denominator:.4f}")
    print(f"  - Expected PGS R^2 = {expected_r2_pgs:.4f}")
    print("\n" + "="*40 + "\n")

    # --- Conclusion ---
    print("Conclusion:")
    print(f"The calculated predictive ability of the PGS (R^2 = {expected_r2_pgs:.4f}) is lower than the total SNP heritability (h2_snp = {h2_snp}).")
    print("This demonstrates that the predictive power of a PGS is constrained by factors like GWAS sample size and is necessarily lower than the theoretical maximum.")

demonstrate_pgs_prediction_limit()