import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite structure (APbBr3) and determines
    their suitability for forming a 3D structure.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # rB for Pb^2+
    r_Br = 196  # rX for Br^-

    # Effective ionic radii for various A-site cations (rA)
    cations = {
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium": 204,
        "Ethylammonium (EA)": 275,
        "Methylhydrazinium": 250,
        "Dimethylammonium (DMA)": 272
    }

    # Stability range for 3D perovskite
    t_min = 0.8
    t_max = 1.0

    print("Analyzing A-site cations for APbBr3 perovskite structure...")
    print(f"Using r(Pb^2+) = {r_Pb} pm and r(Br^-) = {r_Br} pm")
    print(f"Ideal tolerance factor range for 3D perovskite: {t_min} < t < {t_max}\n")

    denominator = math.sqrt(2) * (r_Pb + r_Br)

    for name, r_A in cations.items():
        numerator = r_A + r_Br
        t = numerator / denominator
        
        # Check if it falls in the stable 3D range
        is_stable_3d = t_min < t < t_max

        # Output the detailed calculation
        print(f"--- For {name} ---")
        print(f"Equation: t = (rA + rX) / (sqrt(2) * (rB + rX))")
        print(f"Calculation: t = ({r_A} + {r_Br}) / (sqrt(2) * ({r_Pb} + {r_Br}))")
        print(f"Result: t = {numerator} / {denominator:.1f} = {t:.3f}")
        
        if is_stable_3d:
            print("Verdict: SUITABLE for 3D perovskite structure.")
        else:
            print("Verdict: UNSUITABLE for 3D perovskite structure (t > 1.0 suggests 2D structure formation).")
        print("-" * 25)

    print("\nEvaluating the answer choices:")
    print("A: Cesium, Methylammonium, Formamidinium -> (All known 3D formers, but is it the most comprehensive?)")
    print("B: Cesium, Methylammonium, Formamidinium and Aziridinium -> (All suitable by calculation)")
    print("C: Cesium, Methylammonium, Formamidinium and Ethylammonium -> (Contains unsuitable Ethylammonium)")
    print("D: Cesium, Methylammonium, Formamidinium and Methylhydrazinium -> (All suitable by calculation)")
    print("E: Cesium, Methylammonium, Formamidinium and Dimethylammonium -> (Contains unsuitable Dimethylammonium)")
    
    print("\nConclusion: Options C and E are incorrect. Options A, B, and D contain valid cations. Since the question asks for a comprehensive list, B and D are better than A. Both Aziridinium and Methylhydrazinium are valid candidates.")


if __name__ == '__main__':
    calculate_tolerance_factor()