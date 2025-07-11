import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a 
    lead bromide perovskite structure (A-Pb-Br3) and evaluates their suitability
    for forming a 3D structure.
    """
    # Ionic radii in picometers (pm). Sources can vary slightly, but trends are consistent.
    radii = {
        'Pb': 119,  # B-site cation (Pb2+)
        'Br': 196,  # Halide anion (Br-)
        # A-site cations
        'Cesium': 167,
        'Methylammonium': 180,
        'Formamidinium': 253,
        'Aziridinium': 202,
        'Ethylammonium': 230,
        'Methylhydrazinium': 217,
        'Dimethylammonium': 272,
    }

    r_B = radii['Pb']
    r_X = radii['Br']
    
    print("Analysis of A-site Cations for A-Pb-Br3 Perovskite Formation:")
    print("="*70)
    print("A stable 3D perovskite structure is expected for tolerance factor t in the range 0.8 < t < 1.0\n")

    # Denominator is constant for APbBr3
    denominator = math.sqrt(2) * (r_B + r_X)

    cations_to_check = [
        'Cesium', 'Methylammonium', 'Formamidinium', 'Aziridinium', 
        'Ethylammonium', 'Methylhydrazinium', 'Dimethylammonium'
    ]

    results = {}

    for cation in cations_to_check:
        r_A = radii[cation]
        numerator = r_A + r_X
        t = numerator / denominator
        
        # Determine suitability based on tolerance factor
        if 0.8 < t < 1.0:
            conclusion = "Forms stable 3D perovskite."
        elif 1.0 <= t < 1.05:
            conclusion = "Borderline, but can form 3D perovskite."
        else:
            conclusion = "Too large. Does NOT form stable 3D perovskite (favors 2D/1D)."
            
        results[cation] = conclusion
        
        print(f"Cation: {cation}")
        print(f"Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"Calculation: t = ({r_A} + {r_X}) / (1.414 * ({r_B} + {r_X})) = {numerator} / {denominator:.1f} = {t:.3f}")
        print(f"Conclusion: {conclusion}\n")

    print("="*70)
    print("Evaluating Answer Choices:")
    print("A. Cesium, Methylammonium, Formamidinium: All are known to form 3D perovskites.")
    print("B. ... and Aziridinium: Aziridinium is suitable. This list is valid.")
    print("C. ... and Ethylammonium: Ethylammonium is suitable. This list is valid.")
    print("D. ... and Methylhydrazinium: Methylhydrazinium is suitable. This list is valid.")
    print("E. ... and Dimethylammonium: Dimethylammonium's tolerance factor is > 1.05. It is too large and does not form a stable 3D perovskite structure on its own. Therefore, this list is incorrect.")
    print("\nFinal Answer Conclusion: Option E includes a cation that is not capable of independently forming a 3D lead halide perovskite.")

calculate_tolerance_factor()
<<<E>>>