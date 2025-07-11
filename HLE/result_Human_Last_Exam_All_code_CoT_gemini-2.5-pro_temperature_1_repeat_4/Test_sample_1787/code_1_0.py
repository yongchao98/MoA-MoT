import math

def solve_transit_brightness_drop():
    """
    Calculates the brightness drop of a brown dwarf due to a planetary transit.
    """
    
    # Given ratios from the problem description
    # ratio_1: angular size of Planet 1 / angular size of Planet 2, from the station
    ratio_1 = 0.8
    # ratio_2: angular size of Planet 1 / angular size of the Brown Dwarf, from Planet 2
    ratio_2 = 0.2

    # From the problem statement and orbital mechanics, we deduce the ratio of the orbital distances.
    # Planet 1 is at its pericenter (d1) in a parabolic orbit.
    # Planet 2 is in a circular orbit (d2).
    # The problem implies a specific relationship where d2 is equal to the semi-latus rectum (p)
    # of Planet 1's orbit. For a parabolic orbit, p = 2 * d1.
    # Thus, d2 = 2 * d1, which means the ratio x = d1 / d2 is 0.5.
    x = 0.5

    # From the angular size observations, we can derive the ratio of the radii of Planet 2 and the Brown Dwarf (R₂/R_BD).
    # The derivation is as follows:
    # 1) R₁/d₁ = 0.8 * R₂/d₂  => R₁/R₂ = 0.8 * (d₁/d₂)
    # 2) R₁/(d₂-d₁) = 0.2 * R_BD/d₂ => R₁/R_BD = 0.2 * (d₂-d₁)/d₂ = 0.2 * (1 - d₁/d₂)
    # Dividing R₁/R_BD by R₁/R₂ gives: R₂/R_BD = (0.2/0.8) * (1 - d₁/d₂) / (d₁/d₂)
    # Substituting x = d₁/d₂: R₂/R_BD = 0.25 * (1/x - 1)
    
    R2_div_RBD = (ratio_2 / ratio_1) * (1/x - 1)

    # The brightness drop is determined by the ratio of the planets' disk areas (f).
    f = R2_div_RBD**2

    # The flux ratio during the transit.
    flux_ratio = 1 - f

    # The brightness drop in magnitudes (Δm).
    delta_m = -2.5 * math.log10(flux_ratio)

    # Output the final calculation, showing each number in the equation.
    print("The final calculation for the brightness drop in magnitudes (Δm) is:")
    print(f"Δm = -2.5 * log10(1 - (R₂/R_BD)²) ")
    print(f"Δm = -2.5 * log10(1 - {R2_div_RBD:.2f}²) ")
    print(f"Δm = -2.5 * log10(1 - {f:.4f}) ")
    print(f"Δm = -2.5 * log10({flux_ratio:.4f}) ")
    print(f"Δm = {delta_m:.3f} magnitudes")


solve_transit_brightness_drop()
<<<0.070>>>