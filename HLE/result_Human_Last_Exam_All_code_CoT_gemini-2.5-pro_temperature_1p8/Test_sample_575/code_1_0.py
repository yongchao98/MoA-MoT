def analyze_heritability():
    """
    Analyzes the relationship between heritability and polygenic scores
    based on the problem statement.
    """
    # 1. Setup based on the problem statement
    H_squared = 0.5  # Broad-sense heritability
    # For illustration, let's assume total phenotypic variance (Vp) is 100 units.
    # The exact value doesn't matter, it's the proportions that are important.
    Vp = 100

    # 2. Calculate total genetic variance (Vg)
    # H² = Vg / Vp  =>  Vg = H² * Vp
    Vg = H_squared * Vp

    print("--- Initial Setup Based on Problem ---")
    print(f"Given broad-sense heritability (H²) = Vg / Vp = {H_squared}")
    print(f"Let's assume a total phenotypic variance (Vp) = {Vp} for demonstration.")
    print(f"This implies the total genetic variance (Vg) = H² * Vp = {H_squared} * {Vp} = {Vg}\n")

    print("--- Understanding the Components of Genetic Variance ---")
    # Total genetic variance (Vg) is composed of:
    # Va: Additive genetic variance (the sum of individual allele effects)
    # Vd: Dominance variance (interaction between alleles at the same gene)
    # Vi: Epistatic variance (interaction between alleles at different genes)
    print("Total genetic variance (Vg) has three components:")
    print("Equation: Vg = Va + Vd + Vi")
    print(f"In our example, {Vg} = Va + Vd + Vi")
    print("Since Vd and Vi are variances, their values must be >= 0.")
    print("This means it is necessarily true that: Va <= Vg")
    print(f"So, in our example, Va must be less than or equal to {Vg}\n")

    print("--- Relating Heritability to a Polygenic Score (PGS) ---")
    # A standard Polygenic Score from GWAS is built on a linear model, so it
    # primarily captures additive genetic variance (Va).
    # Narrow-sense heritability (h²) is the proportion of variance from additive effects.
    # The maximum variance a perfect PGS can explain is therefore equal to h².
    print("A Polygenic Score (PGS) predicts a phenotype by summing the effects of many genes.")
    print("Standard GWAS-based PGSs capture the additive genetic effects (Va).")
    print("The proportion of variance due to additive effects is the narrow-sense heritability (h²).")
    print("Equation: h² = Va / Vp")
    print("Therefore, the maximum variance a PGS can explain is h².\n")


    print("--- The Core Conclusion ---")
    # We established that Va <= Vg.
    # If we divide both sides by Vp, the inequality holds: Va / Vp <= Vg / Vp
    # Substituting the definitions of h² and H² gives us the key relationship.
    print("Starting from the necessary truth that Va <= Vg...")
    print("Divide by Vp: (Va / Vp) <= (Vg / Vp)")
    print("Substitute definitions: h² <= H²")
    print(f"Using the number from the problem: h² <= {H_squared}")
    print("Since the variance a PGS can explain is at most h², it logically follows:")
    print(f"Final Equation: Variance_explained_by_PGS <= h² <= {H_squared}")
    print("This shows the PGS can, at most, explain 50% of the variance.\n")

    print("--- Analysis of Answer Choices ---")
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"   - Our final equation (Var_PGS <= {H_squared}) directly supports this.")
    print("   - This statement is NECESSARILY TRUE.\n")

    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("   - This assumes h² = H² = 0.5, which happens only if non-additive effects (Vd + Vi) are zero.")
    print("   - It's possible that Vd or Vi are > 0. For example, if Va=30, Vd=10, Vi=10, then Vg=50 (correct).")
    print(f"   - In that case, h² = Va / Vp = 30 / 100 = 0.3. A perfect PGS would explain 30%, not 50%.")
    print("   - This statement is NOT necessarily true.\n")

    print("C. ...the polygenic score ... will not approach a variance explained of 50% due to non-linear effects...")
    print("   - This assumes that non-additive effects (Vd + Vi) MUST be greater than zero.")
    print("   - However, it's possible that for this trait, Vd=0 and Vi=0.")
    print(f"   - In that case, Va = Vg = {Vg}, and h² = Va / Vp = {Vg} / {Vp} = {H_squared}.")
    print("   - If h²=0.5, a perfect PGS would approach 50%. So this statement is NOT necessarily true.\n")
    
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   - The broad-sense heritability H² = 0.5 already sets the ceiling. It means 50% of variance is genetic and 50% is non-genetic (environmental, epigenetic, etc.).")
    print(f"   - The constraint h² <= {H_squared} is because h² is a component of H², not because of what makes up the non-genetic part.")
    print("   - This statement is NOT necessarily true.\n")

if __name__ == '__main__':
    analyze_heritability()