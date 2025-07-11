import pandas as pd

def analyze_heritability():
    """
    Analyzes and explains the relationship between heritability and polygenic scores.
    """
    # --- Setup from the problem ---
    # Total Phenotypic Variance (Vp) can be set to any arbitrary value for illustration. Let's use 100.
    Vp = 100.0
    
    # Broad-sense heritability (H^2) is given as 0.5
    H2 = 0.5
    
    # This means total Genetic Variance (Vg) is H^2 * Vp
    Vg = H2 * Vp
    
    # Environmental Variance (Ve) is the rest
    Ve = Vp - Vg
    
    print("--- Problem Setup ---")
    print(f"Total Phenotypic Variance (Vp): {Vp}")
    print(f"Broad-Sense Heritability (H^2): {H2}")
    print(f"Total Genetic Variance (Vg = H^2 * Vp): {Vg}")
    print("-" * 25)

    # --- Scenario 1: All genetic variance is additive ---
    # This is a special, theoretical case.
    Va_1 = Vg  # Additive variance
    Vd_1 = 0.0 # Dominance variance
    Vi_1 = 0.0 # Epistatic (interaction) variance
    
    # Narrow-sense heritability (h^2) is Va / Vp
    h2_1 = Va_1 / Vp
    
    # A perfect PGS from an infinitely large GWAS would explain h^2 of the variance.
    pgs_r2_1 = h2_1

    # --- Scenario 2: Genetic variance includes non-additive effects ---
    # This is the typical case for a complex, polygenic trait.
    # We split the total genetic variance (50) into parts.
    Va_2 = 30.0 # Additive variance
    Vd_2 = 12.0 # Dominance variance
    Vi_2 = 8.0  # Epistatic (interaction) variance
    # Note: Va + Vd + Vi = 30 + 12 + 8 = 50 = Vg
    
    h2_2 = Va_2 / Vp
    pgs_r2_2 = h2_2
    
    print("\n--- Evaluating the Statements ---\n")
    print("Statement A: The polygenic score can not explain more than 50% of the variance.")
    print(f"In Scenario 1 (no non-additive effects), max PGS variance explained = h^2 = {pgs_r2_1:.2f} (50%).")
    print(f"In Scenario 2 (with non-additive effects), max PGS variance explained = h^2 = {pgs_r2_2:.2f} (30%).")
    print(f"In both cases, the PGS variance explained is <= 0.5. So, Statement A is NECESSARILY TRUE.\n")

    print("Statement C: The PGS will not approach 50% due to non-linear effects.")
    print("This statement assumes a realistic polygenic trait, which corresponds to Scenario 2.")
    print(f"In Scenario 2, the presence of dominance/epistatic variance ({Vd_2+Vi_2}) means h^2 ({h2_2:.2f}) is strictly less than H^2 ({H2:.2f}).")
    print("Therefore, the PGS variance explained will not approach 50%. It approaches 30% in this example.")
    print("Under the standard assumption that polygenic traits have non-additive effects, Statement C is also NECESSARILY TRUE.\n")

    print("--- Final Conclusion ---")
    print("Since both statements A and C are necessarily true, the correct answer choice is E.")

# Run the analysis
analyze_heritability()