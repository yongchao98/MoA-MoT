def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given facts.
    """
    # Let's assume a total phenotypic variance (Vp) of 100 for easy calculation.
    Vp = 100

    # Given information from the problem
    H2_given = 0.5  # Broad-sense heritability

    # From H2 = Vg / Vp, we can calculate the total genetic variance (Vg)
    Vg = H2_given * Vp
    # Environmental variance (Ve) is the remainder
    Ve = Vp - Vg

    print("--- Problem Setup ---")
    print(f"Given Broad-Sense Heritability (H²) = {H2_given}")
    print(f"Assuming Phenotypic Variance (Vp) = {Vp}")
    print(f"This means Total Genetic Variance (Vg) = {Vg}")
    print(f"And Environmental Variance (Ve) = {Ve}")
    print("\n" + "="*40 + "\n")

    # --- Scenario 1: Purely Additive Genetic Architecture ---
    print("--- Scenario 1: Purely Additive Genetics ---")
    Va_1 = Vg  # All genetic variance is additive
    Vd_1 = 0
    Vi_1 = 0
    h2_1 = Va_1 / Vp
    max_R2_PGS_1 = h2_1
    print(f"Here, we assume all genetic variance is additive (Vd=0, Vi=0).")
    print(f"Additive Variance (Va) = {Va_1}")
    print(f"Narrow-Sense Heritability (h²) = Va / Vp = {Va_1}/{Vp} = {h2_1}")
    print(f"Max variance explained by PGS (R²_PGS) = h² = {max_R2_PGS_1}")
    print("\n")

    # --- Scenario 2: Non-Additive Effects Exist ---
    print("--- Scenario 2: With Non-Additive Genetics ---")
    Va_2 = Vg * 0.6  # Assume 60% of genetic variance is additive
    Vd_2 = Vg * 0.3  # Assume 30% is dominance
    Vi_2 = Vg * 0.1  # Assume 10% is epistatic
    h2_2 = Va_2 / Vp
    max_R2_PGS_2 = h2_2
    print(f"Here, we assume non-additive effects exist (Vd > 0 or Vi > 0).")
    print(f"Additive Variance (Va) = {Va_2}")
    print(f"Dominance Variance (Vd) = {Vd_2}")
    print(f"Epistatic Variance (Vi) = {Vi_2}")
    print(f"Narrow-Sense Heritability (h²) = Va / Vp = {Va_2}/{Vp} = {h2_2}")
    print(f"Max variance explained by PGS (R²_PGS) = h² = {max_R2_PGS_2}")
    print("\n" + "="*40 + "\n")

    # --- Evaluate Each Statement ---
    print("--- Evaluating Statements ---")
    print("A statement is 'necessarily true' only if it holds for ALL valid scenarios.\n")

    # Statement A
    A_scen1 = max_R2_PGS_1 <= 0.5
    A_scen2 = max_R2_PGS_2 <= 0.5
    print("A. The polygenic score can not explain more than 50% of the variance.")
    print(f"   - Scenario 1 (R²_PGS={max_R2_PGS_1}): Is {max_R2_PGS_1} <= 0.5? {A_scen1}")
    print(f"   - Scenario 2 (R²_PGS={max_R2_PGS_2}): Is {max_R2_PGS_2} <= 0.5? {A_scen2}")
    print(f"   - Conclusion: Statement A is NECESSARILY TRUE.\n")

    # Statement B
    B_scen1 = max_R2_PGS_1 == 0.5
    B_scen2 = max_R2_PGS_2 == 0.5
    print("B. The polygenic score will approach a variance explained of 50%.")
    print(f"   - Scenario 1 (R²_PGS={max_R2_PGS_1}): Does {max_R2_PGS_1} approach 0.5? {B_scen1}")
    print(f"   - Scenario 2 (R²_PGS={max_R2_PGS_2}): Does {max_R2_PGS_2} approach 0.5? {B_scen2}")
    print(f"   - Conclusion: Statement B is NOT necessarily true.\n")

    # Statement C
    C_scen1 = max_R2_PGS_1 < 0.5
    C_scen2 = max_R2_PGS_2 < 0.5
    print("C. The PGS will not approach 50% due to non-linear effects.")
    print(f"   - Scenario 1 (R²_PGS={max_R2_PGS_1}): Is {max_R2_PGS_1} not approaching 0.5? {C_scen1}")
    print(f"   - Scenario 2 (R²_PGS={max_R2_PGS_2}): Is {max_R2_PGS_2} not approaching 0.5? {C_scen2}")
    print(f"   - Conclusion: Statement C is NOT necessarily true (fails in Scenario 1).\n")
    
    # Statement D
    # This statement concerns the value of h², independent of epigenetics' contribution to Ve.
    D_scen1 = h2_1 < 0.5
    D_scen2 = h2_2 < 0.5
    print("D. The existence of any epigenetic effects would limit h² to less than 0.5.")
    print("   This is evaluated by checking if h² must be strictly < 0.5.")
    print(f"   - Scenario 1 (h²={h2_1}): Is {h2_1} < 0.5? {D_scen1}")
    print(f"   - Scenario 2 (h²={h2_2}): Is {h2_2} < 0.5? {D_scen2}")
    print(f"   - Conclusion: Statement D is NOT necessarily true (fails in Scenario 1).\n")

# Run the analysis
analyze_heritability_statements()