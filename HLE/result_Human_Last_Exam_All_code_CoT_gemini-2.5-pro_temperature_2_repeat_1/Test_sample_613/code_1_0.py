import math

def calculate_perovskite_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide (APbBr3) perovskite structure and evaluates their suitability.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site cation: Lead (Pb^2+)
    r_Br = 196  # X-site anion: Bromide (Br^-)

    # Effective radii of A-site cations in pm
    # Note: Cesium is inorganic but included as a baseline from the choices.
    cation_radii = {
        "Cesium": 167,
        "Methylammonium": 217,
        "Formamidinium": 253,
        "Aziridinium": 234,
        "Ethylammonium": 274,
        "Methylhydrazinium": 217,
        "Dimethylammonium": 272
    }
    
    print("Calculating Goldschmidt Tolerance Factor (t) for APbBr3 Perovskites...")
    print("Formula: t = (r_A + r_Br) / (sqrt(2) * (r_Pb + r_Br))\n")
    
    # Pre-calculate the denominator as it's constant for all cations here
    denominator = math.sqrt(2) * (r_Pb + r_Br)

    for cation, r_A in cation_radii.items():
        numerator = r_A + r_Br
        t = numerator / denominator
        
        # Determine stability based on the tolerance factor
        if 0.8 < t <= 1.0:
            stability = "Forms stable 3D structure."
        elif t > 1.0 and t < 1.1:
            stability = "Likely forms stable 3D structure (often cubic phase)."
        else:
            stability = "Unlikely to form 3D structure (tends to form 2D)."

        print(f"Cation: {cation} (r_A = {r_A} pm)")
        print(f"t = ({r_A} + {r_Br}) / (sqrt(2) * ({r_Pb} + {r_Br}))")
        print(f"t = {numerator} / {denominator:.2f} = {t:.3f}")
        print(f"Assessment: {stability}\n")

calculate_perovskite_tolerance_factor()