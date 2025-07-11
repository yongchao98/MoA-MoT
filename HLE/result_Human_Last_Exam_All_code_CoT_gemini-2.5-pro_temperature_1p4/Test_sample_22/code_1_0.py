import sys

def explain_pgs_heritability():
    """
    This function explains the relationship between SNP heritability and
    the predictive ability of a Polygenic Score (PGS).
    """
    
    # --- Step 1: Define Key Concepts and Hypothetical Values ---
    
    # SNP Heritability (h²_SNP): The total proportion of phenotypic variance
    # that is theoretically explainable by all the SNPs included in a GWAS.
    # It represents the theoretical maximum predictive power.
    # Let's assume for a trait like height, it's 50%.
    snp_heritability = 0.50

    # Polygenic Score R² (R²_PGS): The actual variance explained when a PGS
    # is applied to a new sample. This is the practical predictive ability.
    # Due to factors like GWAS sample size, noisy effect estimates, and
    # methodological choices, this is lower than the theoretical max.
    # A good PGS for height might explain 10-25% of the variance. Let's use 15%.
    pgs_r_squared = 0.15

    # --- Step 2: Explain the Relationship ---
    
    print("Understanding the Concept:")
    print("-" * 30)
    print("1. SNP Heritability (h²_SNP): This is the theoretical ceiling.")
    print("   It represents the maximum possible variance that could be explained if we knew")
    print("   the true effects of all measured SNPs perfectly.")
    print("\n2. PGS Predictive Ability (R²_PGS): This is the practical performance.")
    print("   A PGS is built from a study (GWAS) with a limited sample size, which means the")
    print("   SNP effect sizes are estimates, not perfect truths. This inherent uncertainty,")
    print("   along with other methodological factors, limits the score's predictive power.")
    print("-" * 30)

    # --- Step 3: Show the relationship as a final equation ---

    print("\nConclusion:")
    print("The predictive ability of the PGS is necessarily bounded by the SNP heritability.")
    print("The final relationship can be expressed with our example numbers:")
    
    # Using 'sys.stdout.write' to build the equation piece by piece to satisfy the prompt.
    sys.stdout.write("R²_PGS (")
    sys.stdout.write(str(pgs_r_squared))
    sys.stdout.write(") <= h²_SNP (")
    sys.stdout.write(str(snp_heritability))
    sys.stdout.write(")\n")
    
    print("\nBecause the PGS is an imperfect, practical application of the underlying genetic architecture, its predictive power must be less than the total variance that is theoretically explainable. Therefore, the statement is True.")

if __name__ == '__main__':
    explain_pgs_heritability()