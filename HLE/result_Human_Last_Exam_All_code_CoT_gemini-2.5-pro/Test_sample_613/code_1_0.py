import math

def calculate_tolerance_factor():
    """
    Calculates and prints the Goldschmidt tolerance factor for various
    A-site cations in lead halide perovskite structures.
    """
    # Ionic Radii in picometers (pm)
    r_B = {'Pb': 119}
    r_X = {'I': 220, 'Br': 196, 'Cl': 181}
    r_A = {
        'Cesium': 188,
        'Methylammonium': 217,
        'Formamidinium': 253,
        'Aziridinium': 178,
        'Ethylammonium': 274,
        'Methylhydrazinium': 217,
        'Dimethylammonium': 272,
    }

    cation_list = [
        'Cesium', 'Methylammonium', 'Formamidinium', 
        'Ethylammonium', 'Dimethylammonium', 'Methylhydrazinium', 'Aziridinium'
    ]
    
    b_cation = 'Pb'
    
    print("--- Goldschmidt Tolerance Factor Calculations for APbX3 Perovskites ---\n")

    for a_cation_name in cation_list:
        ra = r_A[a_cation_name]
        rb = r_B[b_cation]
        
        print(f"--- Evaluating A-site cation: {a_cation_name} (r_A = {ra} pm) ---\n")
        
        for x_anion_name, rx in r_X.items():
            
            numerator = ra + rx
            denominator_part = rb + rx
            denominator = math.sqrt(2) * denominator_part
            t_factor = numerator / denominator

            # Determine stability conclusion
            if 0.8 < t_factor < 1.0:
                conclusion = "Forms stable 3D perovskite (0.8 < t < 1.0)"
            elif t_factor >= 1.0:
                conclusion = "Does NOT form 3D perovskite (t >= 1.0, cation too large, forms 2D/other phases)"
            else: # t_factor <= 0.8
                conclusion = "Does NOT form 3D perovskite (t <= 0.8, cation too small, unstable)"

            # Print the detailed calculation
            print(f"For {a_cation_name} Lead {x_anion_name}ide (A={a_cation_name}, B={b_cation}, X={x_anion_name}):")
            print(f"t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
            print(f"t = ({ra} + {rx}) / (1.414 * ({rb} + {rx}))")
            print(f"t = {numerator} / (1.414 * {denominator_part})")
            print(f"t = {t_factor:.3f}")
            print(f"Conclusion: {conclusion}\n")
        print("-" * 60)

calculate_tolerance_factor()