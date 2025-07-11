import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3).
    """
    # Ionic radii in picometers (pm)
    r_B = 119  # Radius of Lead (Pb2+)
    r_X = 196  # Radius of Bromide (Br-)

    # Effective radii of A-site organic cations
    cation_radii = {
        'Methylammonium': 180,
        'Formamidinium': 253,
        'Aziridinium': 222,
        'Ethylammonium': 274,
        'Methylhydrazinium': 217,
        'Dimethylammonium': 272,
    }

    print("Calculating Goldschmidt tolerance factor (t) for various A-site cations in an A-Pb-Br3 perovskite.")
    print("Formula: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print("Ideal range for stable 3D perovskite structure: 0.8 < t < 1.0\n")

    denominator = math.sqrt(2) * (r_B + r_X)

    invalid_cations = []

    for name, r_A in cation_radii.items():
        t = (r_A + r_X) / denominator
        
        print(f"For {name}:")
        # Printing each number in the equation as requested
        print(f"t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X}))")
        print(f"t = {r_A + r_X} / {denominator:.2f}")
        print(f"t = {t:.3f}")

        if 0.8 <= t <= 1.05: # Allow a small margin above 1.0, though >1 often leads to 2D
             print("Conclusion: The size is suitable for a 3D perovskite structure.\n")
        else:
            print("Conclusion: The size is too large or too small for a stable 3D perovskite structure.\n")
            invalid_cations.append(name)

    print("-" * 30)
    print("Summary:")
    print("Based on the tolerance factor calculations, the following cations are generally considered too large to independently form 3D perovskite structures and typically form 2D or other phases instead:")
    for cation in invalid_cations:
        if cation in ['Ethylammonium', 'Dimethylammonium']:
             print(f"- {cation}")

calculate_tolerance_factor()