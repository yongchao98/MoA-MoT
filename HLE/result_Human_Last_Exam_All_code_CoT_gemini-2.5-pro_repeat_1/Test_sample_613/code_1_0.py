import math

def calculate_tolerance_factor(r_A, r_B, r_X, cation_name):
    """Calculates and prints the Goldschmidt tolerance factor equation and result."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    if denominator == 0:
        t = float('inf')
    else:
        t = numerator / denominator

    suitability = ""
    # The upper limit for 3D perovskites is slightly flexible, but values > 1.05 almost always lead to 2D structures.
    if 1.05 > t > 0.8:
        suitability = "Suitable for 3D structure."
    elif t >= 1.05:
        suitability = "Too large; tends to form 2D or non-perovskite structures."
    else: # t <= 0.8
        suitability = "Too small; tends to form non-perovskite structures."
    
    print(f"Cation: {cation_name}")
    # The final code needs to output each number in the final equation
    print(f"Equation: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X}))")
    print(f"Result: t = {t:.3f}")
    print(f"Assessment: {suitability}\n")

def main():
    """
    Main function to evaluate cations for perovskite structure formation.
    """
    print("Evaluating A-site cations for APbBr3 perovskite formation using the Goldschmidt tolerance factor.")
    print("A stable 3D structure is typically formed when 0.8 < t < 1.0.\n")

    # --- Ionic Radii (in picometers, pm) ---
    # B-site cation: Lead (Pb2+)
    r_B_pb = 119
    # X-site anion: Bromide (Br-)
    r_X_br = 196

    # A-site cations from the answer choices
    # Note: Effective ionic radii for organic cations are based on commonly cited values.
    cations = {
        "Cesium": 167,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 208,
        "Ethylammonium": 274,
        "Methylhydrazinium": 217,
        "Dimethylammonium": 272
    }

    for name, r_A in cations.items():
        calculate_tolerance_factor(r_A, r_B_pb, r_X_br, name)

if __name__ == "__main__":
    main()