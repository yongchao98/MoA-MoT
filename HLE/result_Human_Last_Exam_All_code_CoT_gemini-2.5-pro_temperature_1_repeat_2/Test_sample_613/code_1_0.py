import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide (APbBr3) perovskite structure to determine their
    capability of forming a 3D crystal lattice.
    """

    # Ionic radii in picometers (pm).
    # For organic cations, these are effective ionic radii.
    r_B = 119  # Radius of Pb2+
    r_X = 196  # Radius of Br-

    # Dictionary of A-site cations and their radii from the answer choices
    cations = {
        'Cesium': 167,
        'Methylammonium': 217,
        'Formamidinium': 253,
        'Aziridinium': 202,
        'Ethylammonium': 274,
        'Dimethylammonium': 272,
        'Methylhydrazinium': 230
    }

    print("Calculating Goldschmidt tolerance factor (t) for APbBr3 perovskites.")
    print("A stable 3D structure is generally formed when 0.8 < t < 1.0.\n")
    
    # constant part of the denominator
    sqrt2 = math.sqrt(2)
    denominator_part = r_B + r_X
    
    # We evaluate the cations relevant to the multiple choice options
    cations_to_evaluate = [
        'Cesium', 'Methylammonium', 'Formamidinium', 'Aziridinium', 
        'Ethylammonium', 'Dimethylammonium'
    ]

    for name in cations_to_evaluate:
        r_A = cations[name]
        
        # Calculate tolerance factor
        numerator = r_A + r_X
        denominator = sqrt2 * denominator_part
        t = numerator / denominator

        # Determine if it's likely to form a 3D structure
        if t > 1.05:
            conclusion = "Unlikely to form 3D structure (cation is too large)."
        elif t > 0.8:
            conclusion = "Likely to form 3D structure."
        else:
            conclusion = "Unlikely to form 3D structure (cation is too small)."

        print(f"For {name} (A = {name}):")
        # Print the equation with all numbers filled in
        print(f"t = (r_A + r_X) / (√2 * (r_B + r_X))")
        print(f"t = ({r_A} + {r_X}) / ({sqrt2:.3f} * ({r_B} + {r_X})) = {t:.3f}")
        print(f"Conclusion: {conclusion}\n")

calculate_tolerance_factor()