import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3).
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site cation: Pb^2+
    r_Br = 196  # X-site anion: Br^-

    # Dictionary of effective ionic radii for A-site cations
    cations = {
        "Cesium (Cs)": 167,
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (AZ)": 210,
        "Ethylammonium (EA)": 274,
        "Methylhydrazinium (MHy)": 217,
        "Dimethylammonium (DMA)": 272,
    }
    
    print("Analysis of A-site Cations for A-Pb-Br3 Perovskite Formation")
    print("="*65)
    print(f"Using Goldschmidt Tolerance Factor: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print(f"With r_B(Pb^2+) = {r_Pb} pm and r_X(Br^-) = {r_Br} pm")
    print(f"Stable 3D perovskite range: 0.8 < t < 1.1")
    print("-" * 65)

    denominator = math.sqrt(2) * (r_Pb + r_Br)

    for name, r_A in cations.items():
        numerator = r_A + r_Br
        t = numerator / denominator
        
        # Determine stability based on tolerance factor
        if 0.8 < t < 1.1:
            stability = "Forms stable 3D perovskite"
        else:
            stability = "Does not form stable 3D perovskite (favors 2D)"
            
        print(f"Cation: {name}")
        print(f"  Radius (r_A): {r_A} pm")
        # The final requirement: output each number in the final equation
        print(f"  Calculation: t = ({r_A} + {r_Br}) / (sqrt(2) * ({r_Pb} + {r_Br})) = {numerator} / {denominator:.2f} = {t:.3f}")
        print(f"  Result: {stability}\n")

calculate_tolerance_factor()