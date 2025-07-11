def heritability_analysis():
    """
    Analyzes and explains the concepts of broad-sense (H^2) and
    narrow-sense (h^2) heritability to answer the user's question.
    """
    print("--- Understanding Heritability: Formulas ---")
    print("Phenotypic Variance (Vp) = Genetic Variance (Vg) + Environmental Variance (Ve)")
    print("Genetic Variance (Vg) = Additive (Va) + Dominance (Vd) + Epistatic (Vi)")
    print("Broad-sense Heritability (H^2) = Vg / Vp")
    print("Narrow-sense Heritability (h^2) = Va / Vp\n")

    # --- Step 1: Model the Rabbit Experiment ---
    print("--- Scenario 1: Rabbit Experiment (Purely Additive Genetics) ---")
    print("The problem states H^2 = 0.75 and genetic variance is entirely additive.")
    # Let's use example values that result in H^2 = 0.75
    # If H^2 = Vg / Vp = 0.75, let Vg = 75 and Vp = 100. Then Ve = Vp - Vg = 25.
    # Since it's purely additive, Va = Vg, and Vd = 0, Vi = 0.
    V_a_rabbit = 75
    V_d_rabbit = 0
    V_i_rabbit = 0
    V_e_rabbit = 25

    # Calculate genetic and phenotypic variance
    V_g_rabbit = V_a_rabbit + V_d_rabbit + V_i_rabbit
    V_p_rabbit = V_g_rabbit + V_e_rabbit

    # Calculate heritabilities
    H2_rabbit = V_g_rabbit / V_p_rabbit
    h2_rabbit = V_a_rabbit / V_p_rabbit

    print(f"Given variance components: Va = {V_a_rabbit}, Vd = {V_d_rabbit}, Vi = {V_i_rabbit}, Ve = {V_e_rabbit}")
    print(f"Calculated Genetic Variance (Vg) = {V_a_rabbit} + {V_d_rabbit} + {V_i_rabbit} = {V_g_rabbit}")
    print(f"Calculated Phenotypic Variance (Vp) = {V_g_rabbit} + {V_e_rabbit} = {V_p_rabbit}")
    print(f"Broad-sense Heritability (H^2) = {V_g_rabbit} / {V_p_rabbit} = {H2_rabbit:.2f}")
    print(f"Narrow-sense Heritability (h^2) = {V_a_rabbit} / {V_p_rabbit} = {h2_rabbit:.2f}")
    print("Conclusion: When genetic variance is purely additive, H^2 and h^2 are equal.\n")

    # --- Step 2: Analyze the cause of difference in another species ---
    print("--- Scenario 2: Another Species (with Epistasis - Option C) ---")
    print("Let's see what happens if we introduce epistatic variance (Vi > 0).")
    # We use the same base values but add epistatic variance.
    V_a_species = 75
    V_d_species = 0
    V_i_species = 15  # Introduce epistasis
    V_e_species = 25

    # Recalculate variances and heritabilities
    V_g_species = V_a_species + V_d_species + V_i_species
    V_p_species = V_g_species + V_e_species
    H2_species = V_g_species / V_p_species
    h2_species = V_a_species / V_p_species
    
    print(f"Given variance components: Va = {V_a_species}, Vd = {V_d_species}, Vi = {V_i_species}, Ve = {V_e_species}")
    print(f"Calculated Genetic Variance (Vg) = {V_a_species} + {V_d_species} + {V_i_species} = {V_g_species}")
    print(f"Calculated Phenotypic Variance (Vp) = {V_g_species} + {V_e_species} = {V_p_species}")
    print(f"Broad-sense Heritability (H^2) = {V_g_species} / {V_p_species} = {H2_species:.2f}")
    print(f"Narrow-sense Heritability (h^2) = {V_a_species} / {V_p_species} = {h2_species:.2f}")
    print("Conclusion: The presence of epistatic variance (Vi) causes H^2 to be greater than h^2. This directly addresses why the heritability measures would differ.\n")

    # --- Step 3: Evaluate Option E ---
    print("--- Analysis of Option E ---")
    print("Option E claims that dominance variance (Vd) has 'no impact on h^2'. Let's test this.")
    print("Starting from the purely additive rabbit model (h^2 = 0.75), we add dominance variance (Vd).")
    V_a_test = 75
    V_d_test = 10 # Introduce dominance
    V_i_test = 0
    V_e_test = 25
    V_g_test = V_a_test + V_d_test + V_i_test
    V_p_test = V_g_test + V_e_test
    h2_test = V_a_test / V_p_test
    print(f"If we add Vd = {V_d_test}, the new total phenotypic variance Vp becomes {V_g_test} + {V_e_test} = {V_p_test}.")
    print(f"The new h^2 = {V_a_test} / {V_p_test} = {h2_test:.2f}")
    print(f"The new h^2 ({h2_test:.2f}) is different from the original h^2 (0.75). Therefore, the claim in option E is incorrect.\n")

    print("--- Final Conclusion ---")
    print("A difference between H^2 and h^2 is caused by non-additive genetic variance (Vd or Vi).")
    print("Option C correctly identifies that epistatic interactions (which create Vi) would cause H^2 and h^2 to differ.")
    print("Other options are incorrect because environmental variance (A) and phenotypic variance (B) affect both measures, and dominance variance (E) does impact h^2.")

heritability_analysis()
<<<C>>>