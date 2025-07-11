def analyze_heritability_problem():
    """
    Analyzes the quantitative genetics problem statement by statement.
    """

    print("### Conceptual Framework ###")
    print("Vp: Total Phenotypic Variance")
    print("Vg: Total Genetic Variance")
    print("Va: Additive Genetic Variance (genes contribute independently)")
    print("Vd: Dominance Variance (interaction between alleles at the same locus)")
    print("Vi: Epistatic Variance (interaction between alleles at different loci)")
    print("Ve: Environmental Variance")
    print("\nFundamental Equations:")
    print("1. Vp = Vg + Ve")
    print("2. Vg = Va + Vd + Vi")
    print("-" * 30)

    # Given information
    H2 = 0.5
    print(f"### Analysis based on Given Information (H² = {H2}) ###")
    print("Broad-sense heritability (H²) is the proportion of phenotypic variance from ALL genetic sources.")
    print(f"H² = Vg / Vp = (Va + Vd + Vi) / Vp = {H2}\n")

    print("Narrow-sense heritability (h²) is the proportion from ADDITIVE genetic sources only.")
    print("h² = Va / Vp\n")

    print("A Polygenic Score (PGS) from a standard GWAS is based on an additive model.")
    print("Therefore, the maximum variance a PGS can explain is h².\n")
    print("-" * 30)

    print("### Evaluating the Statements ###\n")

    # --- Statement A ---
    print("▶ Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("The variance explained by a PGS is limited by h². By definition, h² ≤ H².")
    print("Final Equation for Statement A:")
    print(f"Variance_Explained_by_PGS ≤ h² ≤ H² = {H2}")
    print("Since H² is 0.5, the PGS can't explain more than 0.5 (50%) of the variance. This is always true.")
    print("Conclusion: Statement A is NECESSARILY TRUE.\n")

    # --- Statement B ---
    print("▶ Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("A perfect PGS explains h² of the variance. It only approaches 50% (H²) if h² = H².")
    print("This would mean that non-additive variance (Vd + Vi) is zero. The problem does not guarantee this.")
    print("Conclusion: Statement B is NOT necessarily true.\n")

    # --- Statement C ---
    print("▶ Statement C: Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects such as dominance.")
    print("This statement argues that a PGS (measuring h²) will fall short of H² because of non-additive effects (Vd + Vi).")
    print("For a complex 'polygenic' trait, the existence of some non-additive effects is a standard biological assumption.")
    print("Under this assumption, (Vd + Vi) > 0.")
    print("Final Equation for Statement C:")
    print("If (Vd + Vi) > 0, then Va < Vg, which implies h² (Va/Vp) < H² (Vg/Vp).")
    print(f"Therefore, h² < {H2}. The PGS would not reach 50%.")
    print("Conclusion: Assuming the trait is truly polygenic and complex, Statement C is TRUE.\n")

    # --- Statement D ---
    print("▶ Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Epigenetic effects are typically part of Environmental Variance (Ve) in this model.")
    print(f"The condition h² ≤ H² = {H2} is always true by definition, not because of epigenetics.")
    print("The presence of Ve (from any source) is already accounted for in the H² calculation.")
    print("Conclusion: Statement D is NOT necessarily true.\n")
    
    print("-" * 30)
    print("### Final Verdict ###")
    print("Statement A is necessarily true from first principles.")
    print("Statement C is true under the standard biological context of a polygenic trait.")
    print("Therefore, the best choice combines A and C.")

if __name__ == '__main__':
    analyze_heritability_problem()
