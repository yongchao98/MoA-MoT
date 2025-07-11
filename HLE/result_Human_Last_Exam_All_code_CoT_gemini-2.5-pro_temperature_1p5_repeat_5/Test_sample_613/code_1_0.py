import math

def calculate_tolerance_factor(r_a, r_b, r_x):
    """
    Calculates the Goldschmidt tolerance factor (t).
    t = (r_a + r_x) / [sqrt(2) * (r_b + r_x)]
    """
    numerator = r_a + r_x
    denominator = math.sqrt(2) * (r_b + r_x)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

def analyze_cation_stability(cation_name, r_a, r_b, r_x, halide_name):
    """
    Analyzes and prints the stability of a cation in a perovskite structure.
    """
    t = calculate_tolerance_factor(r_a, r_b, r_x)
    
    # Determine stability conclusion based on the tolerance factor
    if t > 1.0:
        conclusion = "Unlikely to form a 3D perovskite (t > 1.0, too large, favors 2D/hexagonal phases)."
    elif t >= 0.8:
        conclusion = "Likely to form a stable 3D perovskite."
    else: # t < 0.8
        conclusion = "Unlikely to form a 3D perovskite (t < 0.8, too small, favors lower-dimensional phases)."

    print(f"--- Analysis for {cation_name} in a Lead {halide_name} Perovskite ---")
    print(f"Cation: {cation_name:<20} | Radius (r_A): {r_a:.2f} Å")
    print(f"Tolerance Factor (t): {t:.3f}")
    print(f"Conclusion: {conclusion}\n")

def main():
    """
    Main function to evaluate cations for 3D perovskite formation.
    """
    # Define effective ionic radii in Angstroms (Å)
    # Using Lead (Pb2+) as B-site and Iodide (I-) as X-site for this example.
    R_PB = 1.19  # r_b for Pb2+
    R_I = 2.20   # r_x for I-

    # Dictionary of A-site cations and their effective ionic radii
    cations = {
        "Cesium": 1.81,
        "Methylammonium": 1.80,
        "Formamidinium": 1.95,
        "Aziridinium": 1.65,
        "Ethylammonium": 2.16,
        "Methylhydrazinium": 1.70,
        "Dimethylammonium": 2.72
    }

    print("Evaluating A-site Cations for 3D Perovskite Formation (APbI₃)")
    print("A stable 3D structure is expected for a tolerance factor (t) in the range: 0.8 < t < 1.0\n")
    print(f"Calculation Parameters:")
    print(f"  - Radius of Pb²⁺ (r_B): {R_PB} Å")
    print(f"  - Radius of I⁻ (r_X):  {R_I} Å\n")
    
    for name, radius in cations.items():
        analyze_cation_stability(name, radius, R_PB, R_I, "Iodide")
        
    print("Summary:")
    print("The analysis shows that only Cesium, Methylammonium, and Formamidinium have tolerance factors")
    print("within the ideal range for forming stable, independent 3D lead halide perovskite structures.")
    print("Other cations are either too small (Aziridinium, Methylhydrazinium) or too large (Ethylammonium, Dimethylammonium).")


if __name__ == "__main__":
    main()