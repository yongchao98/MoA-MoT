import math

def calculate_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite structure (APbBr3).
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site: Pb^2+
    r_Br = 196  # X-site: Br^-

    # Effective ionic radii for potential A-site cations
    a_site_cations = {
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 203,
        "Ethylammonium": 274,
        "Methylhydrazinium": 217,
        "Dimethylammonium": 272
    }
    
    # Constant part of the denominator: sqrt(2) * (r_B + r_X)
    denominator_val = math.sqrt(2) * (r_Pb + r_Br)

    print("Calculating Goldschmidt Tolerance Factor (t) for various A-Pb-Br3 perovskites.")
    print("Formula: t = (r_A + r_Br) / (sqrt(2) * (r_Pb + r_Br))")
    print("-" * 70)
    print(f"Using radii: r_Pb = {r_Pb} pm, r_Br = {r_Br} pm\n")

    # Iterate through cations, calculate and print results
    for name, r_A in a_site_cations.items():
        numerator_val = r_A + r_Br
        t = numerator_val / denominator_val
        
        # Determine stability based on tolerance factor
        if 0.8 < t < 1.05: # A slightly relaxed upper bound is sometimes used
            stability = "Likely to form a stable 3D perovskite."
        else:
            stability = "Unlikely to form a stable 3D perovskite (too large or small)."
        
        print(f"For {name} (r_A = {r_A} pm):")
        # Print the full equation with numbers
        print(f"  t = ({r_A} + {r_Br}) / (1.414 * ({r_Pb} + {r_Br}))")
        print(f"  t = {numerator_val} / {denominator_val:.1f}")
        print(f"  t = {t:.3f}")
        print(f"  Result: {stability}\n")
        
calculate_tolerance_factor()