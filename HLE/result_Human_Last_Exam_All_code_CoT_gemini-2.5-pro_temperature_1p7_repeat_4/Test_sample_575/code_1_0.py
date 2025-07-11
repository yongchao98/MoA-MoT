def analyze_heritability():
    """
    Analyzes the relationship between heritability and polygenic scores
    based on the provided information.
    """

    # --- Step 1: Define knowns and concepts ---
    # H2 is the broad-sense heritability. It's the ratio of total genetic variance (Vg)
    # to total phenotypic variance (Vp).
    H2 = 0.5

    # Vg itself is composed of additive (Va), dominance (Vd), and epistatic (Vi) variance.
    # Vg = Va + Vd + Vi

    # h2 is the narrow-sense heritability, which only considers additive variance.
    # h2 = Va / Vp

    # A standard Polygenic Score (PGS) is an additive linear model. Its predictive power
    # (variance explained) is capped by the narrow-sense heritability, h2.
    # Variance_Explained_PGS <= h2

    print("--- Problem Analysis ---")
    print("Given Information:")
    print(f"1. Broad-sense heritability (H2 = Vg / Vp) = {H2}")
    print("2. A Polygenic Score (PGS) is built from GWAS data (an additive model).")
    print("\nKey Relationships:")
    print("  - Total Genetic Variance (Vg) = Additive (Va) + Dominance (Vd) + Epistatic (Vi) Variance")
    print("  - Max Variance Explained by PGS is equal to Narrow-Sense Heritability (h2 = Va / Vp)")

    # --- Step 2: The Logical Deduction ---
    print("\n--- Logical Deduction ---")
    print("By definition, variance components cannot be negative (e.g., Va >= 0, Vd >= 0).")
    print("Therefore, the additive component (Va) must be less than or equal to the total genetic variance (Vg).")
    print("Equation: Va <= Vg")
    print("\nDividing both sides by the total phenotypic variance (Vp) to get heritability values:")
    print("Equation: Va / Vp <= Vg / Vp")
    print("This translates to: h2 <= H2")

    # --- Step 3: Apply Given Values ---
    print("\n--- Conclusion from Deduction ---")
    print(f"We are given that H2 = {H2}.")
    print(f"Therefore, it must be true that h2 <= {H2}.")
    print("\nSince the maximum variance a PGS can explain is h2, it follows that:")
    final_equation = f"Variance Explained by PGS <= h2 <= {H2}"
    print(f"Final Equation: {final_equation}")
    print(f"This means the polygenic score cannot explain more than {H2 * 100}% of the phenotypic variance.")
    print("\nThis makes statement A necessarily true.")


if __name__ == "__main__":
    analyze_heritability()
<<<A>>>