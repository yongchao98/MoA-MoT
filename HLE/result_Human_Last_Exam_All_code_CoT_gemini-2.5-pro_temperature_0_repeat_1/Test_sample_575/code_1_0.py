def analyze_heritability():
    """
    Analyzes the relationship between heritability and polygenic score (PGS) variance explained.
    """
    # --- Given Information ---
    # For simplicity, let's assume total phenotypic variance is 1.0
    # Heritability is a ratio, so the absolute value doesn't matter.
    Vp = 1.0
    H_squared = 0.5  # Broad-sense heritability

    # --- Derived Quantities ---
    # Total genetic variance (Vg) is H^2 * Vp
    Vg = H_squared * Vp

    print("--- Foundational Concepts ---")
    print(f"Total Phenotypic Variance (Vp) is assumed to be: {Vp}")
    print(f"Broad-Sense Heritability (H^2) is given as: {H_squared}")
    print(f"This means the total Genetic Variance (Vg) is Vp * H^2 = {Vp} * {H_squared} = {Vg}")
    print("\nGenetic Variance (Vg) is composed of Additive (Va), Dominance (Vd), and Epistatic (Vi) variance.")
    print("So, Vg = Va + Vd + Vi")
    print("\nA standard Polygenic Score (PGS) from GWAS primarily captures Additive variance (Va).")
    print("The maximum variance a PGS can explain is the Narrow-Sense Heritability (h^2), where h^2 = Va / Vp.")
    print("-" * 30)

    print("\n--- Evaluating the Core Relationship ---")
    print("By definition, the additive variance (Va) can't be larger than the total genetic variance (Vg).")
    print("Therefore, the following inequality must always be true:")
    print("Va <= Vg")
    print("\nIf we divide both sides by the total phenotypic variance (Vp), we get:")
    print("(Va / Vp) <= (Vg / Vp)")
    print("\nSubstituting the definitions of h^2 and H^2:")
    print("h^2 <= H^2")
    print(f"\nSince the variance explained by a PGS is h^2, and we know H^2 = {H_squared}:")
    print(f"Variance explained by PGS <= {H_squared}")
    print("\nThis means a polygenic score cannot explain more than 50% of the variance in the phenotype.")
    print("This directly supports statement A as being necessarily true.")

    print("-" * 30)
    print("\n--- Why Other Statements Are Not Necessarily True ---")
    print("Statement B says the PGS will approach 50%. This only happens if h^2 = H^2, meaning Vd=0 and Vi=0. We can't assume this.")
    print("Statement C says the PGS will NOT approach 50%. This only happens if Vd>0 or Vi>0. We can't assume this either.")
    print("Statement D is a red herring. The mathematical limit h^2 <= H^2 holds regardless of other biological effects like epigenetics.")

if __name__ == '__main__':
    analyze_heritability()
<<<A>>>