def calculate_heritability(Va, Vd_plus_Vi, Ve, case_name):
    """
    Calculates and prints heritability components based on variance inputs.

    Args:
        Va (float): Additive genetic variance.
        Vd_plus_Vi (float): Non-additive genetic variance (Dominance + Epistasis).
        Ve (float): Environmental variance.
        case_name (str): The name of the scenario being tested.
    """
    Vg = Va + Vd_plus_Vi
    Vp = Vg + Ve
    
    # Handle division by zero, though not expected with the setup
    if Vp == 0:
        print("Total phenotypic variance is zero, cannot calculate heritability.")
        return
        
    H2 = Vg / Vp  # Broad-sense heritability
    h2 = Va / Vp  # Narrow-sense heritability (max variance explained by PGS)

    print(f"--- {case_name} ---")
    print(f"Given Variances:")
    print(f"  Additive Genetic Variance (Va)      = {Va}")
    print(f"  Non-Additive Genetic Variance (Vd+Vi) = {Vd_plus_Vi}")
    print(f"  Environmental Variance (Ve)         = {Ve}\n")
    
    print(f"Calculated Components:")
    print(f"  Total Genetic Variance (Vg = Va + Vd+Vi)        = {Va} + {Vd_plus_Vi} = {Vg}")
    print(f"  Total Phenotypic Variance (Vp = Vg + Ve)        = {Vg} + {Ve} = {Vp}\n")

    print(f"Heritability Calculations:")
    print(f"Broad-Sense Heritability (H^2 = Vg / Vp):")
    print(f"  {Vg} / {Vp} = {H2:.2f}")

    print(f"Narrow-Sense Heritability (h^2 = Va / Vp):")
    print(f"  {Va} / {Vp} = {h2:.2f}")
    print(f"In this case, h^2 {'=' if h2 == H2 else '<'} H^2")
    print("-" * (len(case_name) + 8) + "\n")

# --- SCENARIO 1 ---
# A theoretical case where all genetic variance is purely additive.
# We set H^2 = 0.5. Let Vg=50, Ve=50. Vp=100.
# If all genetic variance is additive, Va=50 and Vd+Vi=0.
calculate_heritability(
    Va=50, 
    Vd_plus_Vi=0, 
    Ve=50, 
    case_name="Case 1: No Non-Additive Genetic Variance"
)

# --- SCENARIO 2 ---
# A more realistic case for a complex trait, where non-additive effects exist.
# We keep H^2 = 0.5. Let Vg=50, Ve=50. Vp=100.
# Now, Vg is split between additive and non-additive parts. Let Va=40, Vd+Vi=10.
calculate_heritability(
    Va=40, 
    Vd_plus_Vi=10, 
    Ve=50, 
    case_name="Case 2: With Non-Additive Genetic Variance"
)
