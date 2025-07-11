def explain_pgs_heritability():
    """
    Explains the relationship between SNP heritability and Polygenic Score (PGS) predictive ability.
    """
    # Define the concepts
    snp_heritability_explanation = "SNP heritability (h²_snp) is the theoretical upper limit of phenotypic variance that can be explained by a set of SNPs."
    pgs_r_squared_explanation = "A Polygenic Score's predictive ability (R²) is the actual variance explained in a new sample using a model built from GWAS data."

    # Explain the relationship
    relationship_explanation = "The PGS R² is necessarily lower than h²_snp because the PGS is built with estimated, not true, SNP effects from a finite-sized GWAS. This introduces statistical error."
    
    # Print the explanation
    print("Step 1: Understanding SNP Heritability")
    print(f"  - {snp_heritability_explanation}")
    print("\nStep 2: Understanding Polygenic Score Predictive Ability")
    print(f"  - {pgs_r_squared_explanation}")
    print("\nStep 3: The Conclusion")
    print(f"  - {relationship_explanation}")
    print("  - Therefore, h²_snp serves as a ceiling that the R² of a PGS can only approach but not exceed.")
    
    # Final Answer
    print("\nThe statement is: True")

if __name__ == "__main__":
    explain_pgs_heritability()