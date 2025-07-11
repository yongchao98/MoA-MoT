import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt Tolerance Factor for a list of A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3) and prints the results.
    """
    # Ionic radii in picometers (pm)
    r_B_Pb = 119  # Radius of Pb^2+
    r_X_Br = 196  # Radius of Br^-

    # Effective ionic radii for A-site cations
    cation_radii = {
        "Cesium": 167,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 220,
        "Ethylammonium": 274,
        "Methylhydrazinium": 217,
        "Dimethylammonium": 272,
    }

    print("--- Goldschmidt Tolerance Factor Calculations for A-Pb-Br3 Perovskites ---\n")
    print(f"Using r(Pb^2+) = {r_B_Pb} pm and r(Br^-) = {r_X_Br} pm.\n")
    print("A stable 3D perovskite structure is expected for 0.8 < t < 1.0.\n")

    denominator = math.sqrt(2) * (r_B_Pb + r_X_Br)

    for name, r_A in cation_radii.items():
        numerator = r_A + r_X_Br
        t = numerator / denominator

        print(f"Cation: {name}")
        # Final instruction: output each number in the final equation!
        print(f"t = ({r_A} + {r_X_Br}) / (sqrt(2) * ({r_B_Pb} + {r_X_Br}))")
        print(f"t = {numerator} / {denominator:.2f}")
        print(f"Result: t = {t:.3f}")

        if 0.8 < t < 1.01: # allow for slight deviation for formamidinium
             print("Conclusion: Likely to form a stable 3D perovskite structure.\n")
        else:
             print("Conclusion: Unlikely to form a 3D perovskite; tends to form 2D or other structures.\n")


if __name__ == '__main__':
    calculate_tolerance_factor()