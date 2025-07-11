import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a
    lead bromide (A-Pb-Br3) perovskite structure to predict structural stability.
    """
    # Ionic radii in picometers (pm).
    # r_B is the B-site cation (Lead)
    # r_X is the X-site anion (Bromide)
    r_B = 119  # Pb^2+
    r_X = 196  # Br^-

    # Effective radii for potential A-site cations (r_A)
    # These are the ions we want to test.
    a_site_cations = {
        "Cesium": 167,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 145,        # Estimated, as it's not a common perovskite cation
        "Ethylammonium": 274,
        "Dimethylammonium": 272,
    }

    print("--- Calculating Goldschmidt Tolerance Factor (t) for APbBr3 Perovskites ---")
    print("Formula: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print("A stable 3D perovskite structure is expected for 0.8 < t < 1.0\n")

    denominator = math.sqrt(2) * (r_B + r_X)

    for name, r_A in a_site_cations.items():
        # Calculate the tolerance factor
        t = (r_A + r_X) / denominator
        
        # Determine stability based on the tolerance factor
        if 0.8 < t < 1.0:
            conclusion = "Likely to form a stable 3D perovskite structure."
        elif t <= 0.8:
            conclusion = "Too small. Favors non-perovskite phases (e.g., yellow phase)."
        else: # t >= 1.0
            conclusion = "Too large. Favors 2D layered structures."

        print(f"Cation: {name}")
        # Print the full equation for transparency
        print(f"  Calculation: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X})) = {t:.2f}")
        print(f"  Conclusion: {conclusion}\n")

# Run the calculation
calculate_tolerance_factor()
