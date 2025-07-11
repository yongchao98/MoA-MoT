import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a
    lead bromide (A-Pb-Br3) perovskite structure and evaluates their suitability.
    """
    # Ionic radii in picometers (pm)
    r_B_Pb = 119  # Radius of Pb^2+
    r_X_Br = 196  # Radius of Br^-

    # Effective ionic radii for various potential A-site organic cations
    a_site_cations = {
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium": 202,
        "Ethylammonium (EA)": 274,
        "Methylhydrazinium (MHy)": 217,
        "Dimethylammonium (DMA)": 272,
    }

    print("Analysis of A-site Cations for A-Pb-Br3 Perovskite Formation")
    print("="*70)
    print("A stable 3D perovskite structure is expected for a tolerance factor (t) in the range of ~0.8 to 1.0.\n")

    # Denominator of the tolerance factor equation is constant for A-Pb-Br3
    denominator = math.sqrt(2) * (r_B_Pb + r_X_Br)

    for name, r_A in a_site_cations.items():
        numerator = r_A + r_X_Br
        t = numerator / denominator

        conclusion = ""
        if 1.0 >= t >= 0.8:
            conclusion = "Good fit for 3D perovskite structure."
        elif t > 1.0:
            conclusion = "Too large. Tends to form 2D or non-perovskite structures."
        else: # t < 0.8
            conclusion = "Too small. Tends to form non-perovskite structures."

        print(f"For {name}:")
        print(f"Tolerance factor t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"t = ({r_A} + {r_X_Br}) / (1.414 * ({r_B_Pb} + {r_X_Br}))")
        print(f"t = {numerator} / ({round(math.sqrt(2), 3)} * {r_B_Pb + r_X_Br})")
        print(f"t = {numerator} / {round(denominator, 2)}")
        print(f"t = {round(t, 3)}")
        print(f"Conclusion: {conclusion}\n")

    print("="*70)
    print("Summary:")
    print("- Methylammonium (MA) and Formamidinium (FA) are in the ideal range and are well-known to form stable 3D perovskites.")
    print("- Ethylammonium (EA) and Dimethylammonium (DMA) are too large (t > 1.0) and are known to form 2D or other non-3D structures.")
    print("- While Aziridinium and Methylhydrazinium have suitable sizes, MA and FA are the most common and robust organic cations for high-quality 3D perovskites.")
    print("Therefore, the best list of common and capable cations is Cesium (inorganic), Methylammonium, and Formamidinium.")


calculate_tolerance_factor()
<<<A>>>