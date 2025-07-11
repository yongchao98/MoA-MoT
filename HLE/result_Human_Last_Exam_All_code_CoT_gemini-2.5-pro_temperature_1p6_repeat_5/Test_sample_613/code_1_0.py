import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor for a perovskite structure."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

def check_stability(t):
    """Checks if the tolerance factor corresponds to a stable 3D perovskite."""
    if 0.8 < t < 1.0:
        return "Likely to form a stable 3D perovskite structure."
    elif t >= 1.0:
        return "Likely to form a 2D or distorted structure (too large)."
    else:
        return "Likely unstable or forms a non-perovskite structure (too small)."

def main():
    # Ionic radii in picometers (pm)
    # Using modern DFT-derived effective radii for organic cations for better accuracy.
    radii = {
        'cations_A': {
            'Cesium': 167,
            'Methylammonium': 217,
            'Formamidinium': 253,
            'Aziridinium': 208,
            'Ethylammonium': 274,
            'Methylhydrazinium': 217,
            'Dimethylammonium': 272
        },
        'cation_B': {'Lead (Pb2+)': 119},
        'anion_X': {'Iodide (I-)': 220}
    }

    r_B = radii['cation_B']['Lead (Pb2+)']
    r_X = radii['anion_X']['Iodide (I-)']
    
    print("Evaluating A-site cations for APbI3 perovskite stability:\n")

    for cation, r_A in radii['cations_A'].items():
        t = calculate_tolerance_factor(r_A, r_B, r_X)
        stability = check_stability(t)
        # We need to print each number in the equation.
        # Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
        print(f"Cation: {cation}")
        print(f"Equation: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X}))")
        print(f"Calculated Tolerance Factor (t): {t:.3f}")
        print(f"Conclusion: {stability}\n")
    
    print("--- Analysis of Answer Choices ---")
    print("Choices C (Ethylammonium) and E (Dimethylammonium) include cations with t > 1.0, making them unlikely to form 3D perovskites independently. They are generally known to form 2D structures.")
    print("Choice A is correct but incomplete, as other organic cations can also form stable 3D perovskites.")
    print("Choices B (Aziridinium) and D (Methylhydrazinium) both include cations with ideal tolerance factors (t ~ 0.9). Both are experimentally verified to form 3D perovskites.")
    print("Since both B and D represent valid extensions to the list, and Methylhydrazinium is a well-studied example with a tolerance factor almost identical to the high-performing Methylammonium, we select it as a representative correct answer.")


if __name__ == '__main__':
    main()
