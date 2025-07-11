import math

def calculate_brightness_drop():
    """
    Solves the astrophysics problem to find the brightness drop in magnitudes.
    
    The solution follows these steps:
    1.  The problem provides two geometric relationships based on angular sizes.
        Let R_BD, R_P1, R_P2 be the radii of the Brown Dwarf, Planet 1, and Planet 2.
        Let r_P1, r_P2 be their orbital distances.
        
        From the space station (at r_SS ≈ R_BD):
        (1) R_P1 / (r_P1 - R_BD) = 0.8 * R_P2 / (r_P2 - R_BD)
        
        From Planet 2:
        (2) R_P1 / (r_P2 - r_P1) = 0.2 * R_BD / r_P2

    2.  Using the assumption that planetary orbits are large compared to the BD's radius
        (r_P1 >> R_BD and r_P2 >> R_BD), equation (1) simplifies to:
        (1') R_P1 / r_P1 = 0.8 * R_P2 / r_P2

    3.  Combining (1') and (2) gives an expression for the ratio of the radii of
        Planet 2 and the Brown Dwarf (R_P2/R_BD) in terms of their orbital radii:
        R_P2 / R_BD = (1/4) * (r_P2/r_P1 - 1)

    4.  The problem states Planet 1 has a parabolic orbit (at pericenter r_P1) and
        Planet 2 has a circular orbit (at r_P2). A physical constraint connecting
        them is to assume they have the same specific angular momentum (h).
        h_P2 = sqrt(G*M*r_P2)  (circular)
        h_P1 = sqrt(2*G*M*r_P1) (parabolic at pericenter)
        Equating them gives: sqrt(G*M*r_P2) = sqrt(2*G*M*r_P1) => r_P2 = 2 * r_P1.

    5.  Substituting r_P2 = 2*r_P1 into the expression from step 3, we can calculate R_P2/R_BD.

    6.  The brightness drop (flux ratio) is (R_P2/R_BD)^2.
    
    7.  The drop in magnitudes is calculated using Δm = -2.5 * log10(1 - flux_ratio).
    """

    # From step 4, the ratio of orbital distances is r_P2 / r_P1 = 2
    orbital_ratio = 2.0
    
    # From step 3, we calculate the ratio of the radii R_P2 / R_BD
    # R_P2/R_BD = (1/4) * (r_P2/r_P1 - 1)
    c1 = 0.8 # Angular size ratio from SS
    c2 = 0.2 # Angular size ratio from P2
    # The factor is (c2 / c1) = 0.2 / 0.8 = 1/4
    radii_ratio_factor = c2 / c1 
    radii_ratio = radii_ratio_factor * (orbital_ratio - 1)
    
    # The brightness drop is the ratio of the areas
    flux_drop_ratio = radii_ratio**2
    
    # The brightness drop in magnitudes is Δm = -2.5 * log10(L_final / L_initial)
    # L_final / L_initial = 1 - flux_drop_ratio
    delta_m = -2.5 * math.log10(1 - flux_drop_ratio)
    
    print("Step 1: Determine the ratio of orbital radii.")
    print(f"Based on the orbital types (parabolic and circular), we deduce r_P2 / r_P1 = {orbital_ratio}")
    print("-" * 30)

    print("Step 2: Calculate the ratio of the planet's radius to the brown dwarf's radius.")
    print(f"The ratio R_P2 / R_BD is calculated as ({c2}/{c1}) * ({orbital_ratio} - 1)")
    print(f"R_P2 / R_BD = {radii_ratio:.4f}")
    print("-" * 30)
    
    print("Step 3: Calculate the brightness drop and convert to magnitudes.")
    print("The brightness drop is the ratio of the areas of the two disks.")
    print(f"Brightness Drop Ratio = (R_P2 / R_BD)^2 = ({radii_ratio:.4f})^2 = {flux_drop_ratio:.4f}")
    
    print("\nThe final equation for the brightness drop in magnitudes is:")
    # The prompt asks to output each number in the final equation.
    print(f"Δm = -2.5 * log10(1 - {flux_drop_ratio:.4f})")
    print(f"Δm = -2.5 * log10({1 - flux_drop_ratio:.4f})")
    
    print("\nFinal Answer:")
    # Result accurate to three decimal places.
    print(f"The brightness drop is {delta_m:.3f} bolometric magnitudes.")

calculate_brightness_drop()