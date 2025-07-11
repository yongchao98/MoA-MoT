def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given information.
    """

    # --- Step 1 & 2: Define Terms and Relationships ---
    # H2 = Broad-sense heritability = Vg / Vp (Total Genetic Variance / Phenotypic Variance)
    # h2 = Narrow-sense heritability = Va / Vp (Additive Genetic Variance / Phenotypic Variance)
    # Vg = Va + Vd + Vi (Genetic = Additive + Dominance + Epistatic Variance)
    # A standard Polygenic Score (PGS) from GWAS is a linear model, so its predictive power
    # (R2_PGS) is theoretically limited by narrow-sense heritability (h2).
    # Therefore, we have the fundamental inequality: R2_PGS <= h2 <= H2.

    # --- Step 3: Use Given Information ---
    H2 = 0.5

    print("--- Analysis of the Problem ---")
    print(f"We are given that broad-sense heritability (H2) is {H2}.")
    print("This means total genetic factors account for 50% of the variance in the phenotype.")
    print("The performance of a polygenic score (R2_PGS) is limited by narrow-sense heritability (h2), which in turn is limited by H2.")
    print("The core relationship is: R2_PGS <= h2 <= H2.\n")

    # --- Step 4: Evaluate Each Statement ---

    # Analysis of Statement A
    print("--- Evaluating Statement A ---")
    print("Statement: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"This means R2_PGS <= {H2}.")
    print(f"From our core relationship, R2_PGS <= h2 <= H2. Since H2 = {H2}, it is always true that R2_PGS <= {H2}.")
    print("This statement sets the absolute maximum ceiling for the PGS, which is the total genetic variance.")
    print("Conclusion: Statement A is necessarily TRUE.\n")

    # Analysis of Statement C
    print("--- Evaluating Statement C ---")
    print("Statement: Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects such as dominance.")
    print("A linear PGS's performance is limited by h2 (additive variance), not H2 (total genetic variance).")
    print("The difference between H2 and h2 is the non-additive (non-linear) variance from dominance and epistasis (gene-gene interactions).")
    print("For a 'polygenic' trait, it is a standard and biologically realistic assumption that such non-linear effects exist, meaning h2 < H2.")
    print(f"If h2 < {H2}, then the linear PGS can never reach the ceiling of {H2}, which is 50%. The reason is precisely the non-linear effects mentioned.")
    print("Conclusion: Statement C correctly describes the limitations of a linear PGS and is considered TRUE in this context.\n")

    # Analysis of Statements B and D
    print("--- Evaluating Statements B and D ---")
    print("Statement B claims the PGS *will* approach 50%, which is contradicted by our analysis for C.")
    print("Statement D claims epigenetics *limits* h2 to less than 0.5. The fact that H2=0.5 already limits h2 to a maximum of 0.5. Epigenetic effects contribute to environmental variance and do not necessarily force h2 to be strictly less than H2.")
    print("Conclusion: Statements B and D are NOT necessarily true.\n")
    
    # --- Step 5: Synthesize ---
    print("--- Final Conclusion ---")
    print("Both statements A and C are correct.")
    print("A is true from mathematical definition. C is true from a standard, practical understanding of polygenic architecture and the limitations of linear models.")
    print("Therefore, the correct choice is E, which states that only choices A and C are correct.")

# Run the analysis
analyze_heritability_statements()
