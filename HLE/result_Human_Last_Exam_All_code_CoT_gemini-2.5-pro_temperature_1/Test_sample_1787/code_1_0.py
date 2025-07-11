import math

def solve_transit_brightness_drop():
    """
    Calculates the brightness drop of a brown dwarf during a transit of Planet 2.
    """

    # Step 1 & 2: Define constants and set up equations from the problem statement.
    # Angular size ratio of P1 to P2 from the space station
    k1 = 0.80
    # Angular size ratio of P1 to the Brown Dwarf from Planet 2
    k2 = 0.20

    # The geometric relations are:
    # 1) R1 / (r1 - R_BD) = k1 * R2 / (r2 - R_BD)
    # 2) R1 / (r2 - r1) = k2 * R_BD / r2
    # Using the approximation r1, r2 >> R_BD, equation (1) simplifies to:
    # 1') R1 / r1 = k1 * R2 / r2
    #
    # From (1'), we get R1 = k1 * R2 * (r1 / r2).
    # Substituting this into (2), we get:
    # (k1 * R2 * r1 / r2) / (r2 - r1) = k2 * R_BD / r2
    # Solving for R2/R_BD gives:
    # R2 / R_BD = (k2 / k1) * (r2 / r1 - 1)

    print("Deriving the ratio of the radii of Planet 2 and the Brown Dwarf (R2 / R_BD)...")
    print(f"Ratio of angular sizes k1 (P1/P2 from SS) = {k1}")
    print(f"Ratio of angular sizes k2 (P1/BD from P2) = {k2}")

    # Step 3: Apply the key physical constraint from orbital mechanics.
    # For a stable system containing a circular orbit (r2) and a parabolic orbit (pericenter r1),
    # the condition r2 >= 2*r1 must be met. We assume the limiting case.
    r_ratio = 2.0  # r2 / r1
    print(f"From orbital stability, the ratio of orbital distances r2/r1 = {r_ratio}")

    # Step 4: Solve for the radius ratio R2 / R_BD.
    ratio_R2_RBD = (k2 / k1) * (r_ratio - 1)
    
    print("\nCalculating the final equation for the radius ratio:")
    print(f"R2 / R_BD = (k2 / k1) * (r2/r1 - 1)")
    print(f"R2 / R_BD = ({k2} / {k1}) * ({r_ratio} - 1)")
    print(f"R2 / R_BD = {ratio_R2_RBD}")

    # Step 5: Calculate the brightness drop in magnitudes.
    # The fractional drop in flux is the ratio of the objects' areas.
    flux_drop_fraction = ratio_R2_RBD**2
    
    # The ratio of flux during transit to the flux outside of transit.
    flux_ratio = 1 - flux_drop_fraction
    
    # The change in magnitude is given by the formula: Δm = -2.5 * log10(flux_ratio).
    delta_m = -2.5 * math.log10(flux_ratio)

    print("\nCalculating the brightness drop in magnitudes (Δm)...")
    print(f"Fractional flux drop = (R2 / R_BD)^2 = ({ratio_R2_RBD})^2 = {flux_drop_fraction}")
    print(f"Flux ratio (F_transit / F_normal) = 1 - {flux_drop_fraction} = {flux_ratio}")
    print("Δm = -2.5 * log10(Flux ratio)")
    print(f"Δm = -2.5 * log10({flux_ratio})")
    print(f"The final brightness drop is {delta_m:.3f} magnitudes.")
    
    return delta_m

if __name__ == '__main__':
    final_answer = solve_transit_brightness_drop()
    # The final answer is requested in a specific format at the end.
    # print(f"<<<{final_answer:.3f}>>>")