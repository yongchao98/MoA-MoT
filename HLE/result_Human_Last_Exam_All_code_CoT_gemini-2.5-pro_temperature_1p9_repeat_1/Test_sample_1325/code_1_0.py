def calculate_heritability(va, vd, ve, species_name):
    """
    Calculates and prints broad-sense (H^2) and narrow-sense (h^2) heritability.
    
    Args:
        va (float): Additive genetic variance.
        vd (float): Dominance genetic variance.
        ve (float): Environmental variance.
        species_name (str): Name of the species/population for printing.
    """
    # Phenotypic variance is the sum of all variance components.
    # Assuming epistatic variance (VI) is 0 for simplicity.
    vi = 0
    vp = va + vd + vi + ve
    
    # Genetic variance is the sum of genetic components.
    vg = va + vd + vi

    # Broad-sense heritability (H^2)
    h_squared_broad = vg / vp
    
    # Narrow-sense heritability (h^2)
    h_squared_narrow = va / vp

    print(f"--- {species_name} ---")
    print(f"Variance Components: VA={va}, VD={vd}, VE={ve}")
    print(f"Total Phenotypic Variance (VP) = {va} + {vd} + {ve} = {vp:.2f}")
    
    # Print the equation for Broad-Sense Heritability (H^2)
    print(f"\nBroad-Sense Heritability (H^2) = (VA + VD) / VP")
    print(f"H^2 = ({va} + {vd}) / {vp:.2f} = {h_squared_broad:.3f}")
    
    # Print the equation for Narrow-Sense Heritability (h^2)
    print(f"\nNarrow-Sense Heritability (h^2) = VA / VP")
    print(f"h^2 = {va} / {vp:.2f} = {h_squared_narrow:.3f}\n")

# Scenario 1: Rabbit population
# Given H^2 = 0.75 and genetic variance is entirely additive (VD = 0).
# Let's set VA = 3. Then H^2 = VA / (VA + VE) => 0.75 = 3 / (3 + VE) => VE = 1.
va_rabbit = 3.0
vd_rabbit = 0.0
ve_rabbit = 1.0
calculate_heritability(va_rabbit, vd_rabbit, ve_rabbit, "Rabbit Population (purely additive)")

# Scenario 2: Another Species with Dominance Variance
# Let's introduce some dominance variance (VD) and see its effect.
# We'll keep VA and VE the same to isolate the effect of VD.
va_other = 3.0
vd_other = 2.0  # Introduce dominance variance
ve_other = 1.0
calculate_heritability(va_other, vd_other, ve_other, "Other Species (with dominance)")
