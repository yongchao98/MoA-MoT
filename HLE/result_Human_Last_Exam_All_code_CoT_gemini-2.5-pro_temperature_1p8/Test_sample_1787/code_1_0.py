import math

def calculate_brightness_drop():
    """
    Calculates the brightness drop of a brown dwarf during a transit of Planet 2.

    The solution follows these steps:
    1. Define the given ratios of angular sizes.
    2. Establish equations relating the radii and orbital distances of the planets and the brown dwarf.
    3. Use the physics of the specified orbits (parabolic and circular) to find the ratio of their orbital distances. A key physical assumption is that their specific angular momenta are equal in this special configuration.
    4. Solve for the ratio of the radius of Planet 2 to the radius of the brown dwarf.
    5. Calculate the fractional drop in brightness.
    6. Convert the brightness drop to a change in magnitude.
    """

    # --- Step 1 & 2: Set up geometric relations ---
    # From observation on the station: θ_P1_S = 0.8 * θ_P2_S
    # This leads to: R_P1 / r_P1 = 0.8 * R_P2 / r_P2
    # Let x = r_P1 / r_P2. Then R_P1 / R_P2 = 0.8 * x
    k1 = 0.8

    # From observation on Planet 2: θ_P1_P2 = 0.2 * θ_BD_P2
    # This leads to: R_P1 / (r_P2 - r_P1) = 0.2 * R_BD / r_P2
    # Rearranging: R_P1 / R_BD = 0.2 * (1 - r_P1 / r_P2) = 0.2 * (1 - x)
    k2 = 0.2

    # --- Step 3: Use orbital mechanics to find x = r_P1 / r_P2 ---
    # We assume the specific angular momenta l=r*v are equal.
    # l_P1 = sqrt(2*G*M*r_P1) for a parabolic orbit at pericenter.
    # l_P2 = sqrt(G*M*r_P2) for a circular orbit.
    # Setting l_P1 = l_P2 => sqrt(2*r_P1) = sqrt(r_P2) => 2*r_P1 = r_P2
    # Therefore, x = r_P1 / r_P2 = 0.5
    x = 0.5
    
    # --- Step 4: Solve for the radii ratio R_P2 / R_BD ---
    # We have two expressions for the ratio R_P1 / R_BD:
    # 1) From R_P1/R_P2 = k1 * x, we get R_P1 = (k1 * x) * R_P2.
    #    So, R_P1 / R_BD = (k1 * x) * (R_P2 / R_BD)
    # 2) We also have R_P1 / R_BD = k2 * (1 - x)
    #
    # Equating them: (k1 * x) * (R_P2 / R_BD) = k2 * (1 - x)
    # R_P2 / R_BD = (k2 * (1 - x)) / (k1 * x)
    
    ratio_r_p2_r_bd = (k2 * (1 - x)) / (k1 * x)
    
    # --- Step 5: Calculate the brightness drop ratio ---
    # Brightness drop is the ratio of the areas.
    brightness_drop_ratio = ratio_r_p2_r_bd ** 2

    # --- Step 6: Calculate the magnitude change ---
    # Δm = -2.5 * log10(1 - brightness_drop_ratio)
    mag_drop = -2.5 * math.log10(1 - brightness_drop_ratio)

    print(f"Let x = r_P1/r_P2.")
    print(f"From the observation at the station, we have R_P1/R_P2 = {k1} * x.")
    print(f"From the observation at Planet 2, we have R_P1/R_BD = {k2} * (1 - x).")
    print(f"From the orbital physics (assuming l1=l2), we find x = {x:.1f}.")
    print("\nWe solve for the ratio of Planet 2's radius to the Brown Dwarf's radius (R_P2/R_BD):")
    print(f"R_P2/R_BD = ({k2} * (1 - {x:.1f})) / ({k1} * {x:.1f}) = {ratio_r_p2_r_bd:.4f}")

    print("\nThe fractional brightness drop is the square of this ratio:")
    print(f"ΔF/F = (R_P2/R_BD)^2 = {ratio_r_p2_r_bd:.4f}^2 = {brightness_drop_ratio:.4f}")
    
    print("\nFinally, the brightness drop in magnitudes (Δm) is calculated as:")
    print(f"Δm = -2.5 * log10(1 - ΔF/F)")
    print(f"Δm = -2.5 * log10(1 - {brightness_drop_ratio:.4f})")
    print(f"Δm = -2.5 * log10({1 - brightness_drop_ratio:.4f})")
    print(f"Δm ≈ {mag_drop:.3f}")
    
    # Final answer format
    print(f"\n<<<The brightness drop is {mag_drop:.3f} bolometric magnitudes.>>>")


calculate_brightness_drop()