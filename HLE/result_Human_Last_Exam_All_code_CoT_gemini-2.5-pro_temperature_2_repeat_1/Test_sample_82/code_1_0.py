def analyze_heritability():
    """
    Analyzes the statements about heritability and polygenic scores.
    """
    # --- Step 1: Define concepts and given information ---
    print("--- Defining Key Concepts and Given Information ---")
    H2 = 0.5
    print(f"Given Broad-Sense Heritability (H_squared) = Vg / Vp = {H2}")
    print("Where: ")
    print("Vp = Total Phenotypic Variance")
    print("Vg = Total Genetic Variance")
    print("Va = Additive Genetic Variance")
    print("Vd = Dominance Variance")
    print("Vi = Epistatic (Interaction) Variance")
    print("Vg = Va + Vd + Vi\n")
    print("Narrow-Sense Heritability (h_squared) = Va / Vp")
    print("A standard Polygenic Score (PGS) from GWAS explains a proportion of variance that approaches h_squared.\n")

    # --- Step 2: Create a numerical example ---
    print("--- Creating a Numerical Example ---")
    Vp = 100.0
    Vg = H2 * Vp
    Ve = Vp - Vg
    print(f"Let's assume the total phenotypic variance (Vp) is {Vp}.")
    print(f"Then total genetic variance (Vg) = {H2} * {Vp} = {Vg}")
    print(f"And environmental variance (Ve) = {Vp} - {Vg} = {Ve}\n")

    # For a complex polygenic trait, we assume non-additive variance exists.
    # Let's assign some plausible values.
    Va = 35.0
    Vd = 10.0
    Vi = 5.0
    print("For a polygenic trait, we assume non-additive effects exist (Vd > 0 or Vi > 0).")
    print(f"Let's assume a breakdown of Vg: Va = {Va}, Vd = {Vd}, Vi = {Vi}")
    print(f"Check: {Va} + {Vd} + {Vi} = {Va + Vd + Vi}, which correctly equals our Vg of {Vg}.\n")

    # --- Step 3: Calculate narrow-sense heritability (h2) ---
    h2 = Va / Vp
    print("--- Calculating Narrow-Sense Heritability (h2) ---")
    print(f"The maximum variance an additive PGS can explain is h2 = Va / Vp")
    print(f"h2 = {Va} / {Vp} = {h2}\n")


    # --- Step 4: Evaluate the statements ---
    print("--- Evaluating the Statements ---\n")

    # Statement A
    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    max_pgs_r2 = h2
    print(f"The PGS can explain at most {max_pgs_r2 * 100:.1f}% of the variance in this example.")
    print(f"In general, PGS variance explained (h2) must be less than or equal to total genetic variance (H2).")
    print(f"Since H2 = {H2}, h2 must be <= {H2}. Thus, a PGS cannot explain more than 50%.")
    print("Result: Statement A is necessarily TRUE.\n")


    # Statement C
    print("Statement C: Given an arbitrarily large GWAS, the polygenic score... will not approach a variance explained of 50% due to gene-gene interactions and... dominance.")
    print(f"In our example, the presence of dominance (Vd={Vd}) and epistasis (Vi={Vi}) makes Va < Vg.")
    print(f"This results in h2 ({h2}) being less than H2 ({H2}).")
    print("Therefore, the PGS approaches a value less than 50%.")
    print("Result: Statement C is considered TRUE under the standard interpretation of a 'polygenic trait'.\n")
    
    # --- Final Conclusion ---
    print("--- Final Conclusion ---")
    print("Statements A and C are the correct conclusions based on the definitions.")
    print("Therefore, the correct choice is E, which states that only choices A and C are correct.")


if __name__ == '__main__':
    analyze_heritability()
<<<E>>>