import math

def calculate_perovskite_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a
    lead bromide (A-Pb-Br3) perovskite structure and evaluates their suitability.
    """
    # Ionic radii in picometers (pm)
    # B-site (Lead) and X-site (Bromide) are constant for this problem.
    r_B_val = 119  # Ionic radius of Pb^2+
    r_X_val = 196  # Ionic radius of Br^-

    # A-site cations and their effective ionic radii
    cations = {
        'Cesium': 167, # Inorganic reference
        'Methylammonium': 180,
        'Formamidinium': 253,
        'Aziridinium': 202,
        'Ethylammonium': 270,
        'Methylhydrazinium': 217,
        'Dimethylammonium': 272
    }

    print("Evaluating A-site cations for 3D Lead Bromide Perovskite (A-Pb-Br3) formation:\n")

    # Define the ions from the answer choices to be checked
    cations_to_check = [
        'Cesium',
        'Methylammonium',
        'Formamidinium',
        'Aziridinium',        # From choice B
        'Ethylammonium',     # From choice C
        'Methylhydrazinium', # From choice D
        'Dimethylammonium'   # From choice E
    ]
    
    # Calculate √2 once
    sqrt_2 = math.sqrt(2)

    for name in cations_to_check:
        r_A_val = cations[name]
        
        # Calculate components of the formula
        numerator = r_A_val + r_X_val
        denominator_B_X = r_B_val + r_X_val
        denominator = sqrt_2 * denominator_B_X
        
        # Calculate the tolerance factor
        t = numerator / denominator
        
        # Determine the verdict based on the tolerance factor range [0.8, 1.0]
        if 0.8 <= t <= 1.0:
            verdict = "Likely to form a 3D perovskite."
        else:
            verdict = "Unlikely to form a 3D perovskite (t is too large/small)."
        
        # Print the detailed equation and result
        print(f"Cation: {name}")
        print(f"Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"Calculation: t = ({r_A_val} + {r_X_val}) / (1.414 * ({r_B_val} + {r_X_val})) = {numerator} / (1.414 * {denominator_B_X}) = {t:.3f}")
        print(f"Verdict: {verdict}\n")

# Run the calculation
calculate_perovskite_tolerance_factor()