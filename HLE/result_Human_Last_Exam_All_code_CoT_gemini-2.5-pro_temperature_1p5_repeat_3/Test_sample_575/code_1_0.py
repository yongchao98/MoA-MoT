def analyze_heritability_problem():
    """
    Analyzes the relationship between heritability and polygenic scores
    to determine the correct statement.
    """

    # Step 1 & 2: Define terms and establish relationships
    H2_given = 0.5

    print("--- Analysis of Heritability and Polygenic Scores ---")
    print("\nStep 1: Understanding the concepts")
    print("Broad-Sense Heritability (H^2): The proportion of total phenotypic variance (Vp) that is due to all genetic variance (Vg).")
    print("  H^2 = Vg / Vp")
    print("\nGenetic Variance (Vg) is composed of:")
    print("  - Additive variance (Va): Due to the average linear effects of alleles.")
    print("  - Dominance variance (Vd): Due to interactions between alleles at the same locus.")
    print("  - Epistatic variance (Vi): Due to interactions between alleles at different loci.")
    print("  So, Vg = Va + Vd + Vi")
    print("\nNarrow-Sense Heritability (h^2): The proportion of phenotypic variance (Vp) due to only ADDITIVE genetic variance (Va).")
    print("  h^2 = Va / Vp")
    print("\nPolygenic Score (PGS): A standard PGS is constructed by summing the linear effects of genetic variants from a GWAS. Therefore, the maximum phenotypic variance it can explain is equal to the narrow-sense heritability (h^2).")

    # Step 3: Establish the key inequality
    print("\nStep 2: Relating H^2 and h^2")
    print("Because Vg = Va + Vd + Vi, and variances cannot be negative, it is always true that Va <= Vg.")
    print("If we divide all terms by Vp, the inequality holds: (Va / Vp) <= (Vg / Vp)")
    print("Substituting the definitions of heritability, we get: h^2 <= H^2")

    # Step 4: Apply the given information from the problem
    print(f"\nStep 3: Applying the given information")
    print(f"The problem states that broad-sense heritability H^2 = {H2_given}")
    print(f"Therefore, we can conclude that h^2 <= {H2_given}.")
    print(f"Since a PGS can explain at most h^2 of the variance, the PGS can explain at most {H2_given*100}% of the variance.")

    # Step 5: Evaluate the answer choices
    print("\nStep 4: Evaluating the answer choices")
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("   -> Our analysis shows PGS_variance <= h^2 <= 0.5. This statement is NECESSARILY TRUE.\n")

    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("   -> This would only be true if h^2 = H^2, which requires non-additive variance (Vd + Vi) to be zero. This is not guaranteed.\n")

    print("C. ...the polygenic score...will not approach...50% due to gene-gene interactions...")
    print("   -> This would only be true if h^2 < H^2, which requires non-additive variance to be greater than zero. This is also not guaranteed.\n")

    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   -> Epigenetic effects contribute to non-genetic variance, which is already accounted for in H^2 = 0.5. This statement is irrelevant to the h^2 <= H^2 relationship.\n")
    
    print("--- Final Conclusion ---")
    print("The only statement that must be true based on the provided information is that the variance explained by the polygenic score cannot exceed the broad-sense heritability.")
    final_equation_explanation = f"Variance_Explained_by_PGS <= h^2 <= H^2 = {H2_given}"
    print(f"Final Equation: {final_equation_explanation}")


if __name__ == "__main__":
    analyze_heritability_problem()