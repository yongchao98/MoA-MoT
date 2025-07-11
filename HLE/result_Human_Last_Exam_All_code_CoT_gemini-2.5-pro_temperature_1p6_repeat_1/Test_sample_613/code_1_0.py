import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for potential A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3).
    """

    # Ionic radii in picometers (pm). Effective radii are used for organic cations.
    radii = {
        'B_cation': {'Pb': 119},
        'X_anion': {'Br': 196},
        'A_cation': {
            # Core cations from the options
            'Methylammonium': 217,
            'Formamidinium': 253,
            # Other cations from the options for evaluation
            'Aziridinium': 236,
            'Ethylammonium': 274,
            'Methylhydrazinium': 230,
            'Dimethylammonium': 272,
            # Inorganic reference cation
            'Cesium': 167
        }
    }

    r_B = radii['B_cation']['Pb']
    r_X = radii['X_anion']['Br']
    sqrt_2 = math.sqrt(2)

    print("Evaluating A-site cations for 3D perovskite formation (A-Pb-Br3):")
    print("Stability Range for 3D Perovskite: 0.8 < t < 1.0\n")

    cations_to_test = [
        'Cesium', 'Methylammonium', 'Formamidinium', 'Ethylammonium',
        'Dimethylammonium', 'Aziridinium', 'Methylhydrazinium'
    ]

    for name in cations_to_test:
        r_A = radii['A_cation'][name]
        
        # Calculate the tolerance factor
        numerator = r_A + r_X
        denominator = sqrt_2 * (r_B + r_X)
        t = numerator / denominator

        print(f"--- For {name} in an A-Pb-Br3 structure ---")
        print(f"Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        # Outputting the numbers in the final equation as requested
        print(f"Numbers:  t = ({r_A} + {r_X}) / ({sqrt_2:.3f} * ({r_B} + {r_X}))")
        print(f"Calculation: t = {numerator} / ({sqrt_2:.3f} * {r_B + r_X}) = {numerator} / {denominator:.3f}")
        print(f"Result:   t = {t:.3f}")

        # Evaluate and conclude
        if t > 1.0:
            print("Conclusion: Cation is too large. Tends to form 2D or non-perovskite structures.\n")
        elif t < 0.8:
            print("Conclusion: Cation is too small. Structure is unstable.\n")
        else:
            print("Conclusion: Cation size is suitable for forming a 3D perovskite structure.\n")


if __name__ == '__main__':
    calculate_tolerance_factor()