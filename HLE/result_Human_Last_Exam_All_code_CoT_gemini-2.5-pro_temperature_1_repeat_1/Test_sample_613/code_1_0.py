import math

def analyze_perovskite_cations():
    """
    Analyzes the suitability of various A-site cations for forming 3D lead halide perovskites
    based on the Goldschmidt tolerance factor.
    """

    # Ionic radii in picometers (pm)
    r_B = {'Pb': 119}  # B-site cation: Lead (Pb2+)
    r_X = {'Br': 196}  # X-site anion: Bromide (Br-) as per the question's example

    # Effective ionic radii for A-site cations
    r_A = {
        'Methylammonium (MA+)': 180,
        'Formamidinium (FA+)': 253,
        'Aziridinium (AZ+)': 202,
        'Ethylammonium (EA+)': 274,
        'Methylhydrazinium (MHy+)': 217,
        'Dimethylammonium (DMA+)': 272
    }

    print("Analyzing A-site Cations for 3D Perovskite Formation (A-Pb-Br3)\n")
    print("The Goldschmidt tolerance factor (t) predicts structural stability.")
    print("A stable 3D structure is typically formed when 0.8 < t < 1.0.\n")
    print(f"{'Cation':<25} {'Prediction'}")
    print("-" * 60)

    # Use Bromide as the representative halide
    halide_name = 'Br'
    r_b_val = r_B['Pb']
    r_x_val = r_X[halide_name]

    cations_to_check = [
        'Methylammonium (MA+)',
        'Formamidinium (FA+)',
        'Ethylammonium (EA+)',
        'Dimethylammonium (DMA+)',
        'Aziridinium (AZ+)',
        'Methylhydrazinium (MHy+)'
    ]

    for cation_name in sorted(cations_to_check):
        r_a_val = r_A[cation_name]
        
        # Goldschmidt tolerance factor equation
        numerator = r_a_val + r_x_val
        denominator = math.sqrt(2) * (r_b_val + r_x_val)
        t_factor = numerator / denominator

        # Prediction based on the calculated t-factor
        # Cations with t > ~1.0 tend to form 2D layered structures.
        if 0.8 < t_factor <= 1.05:  # A relaxed upper bound, as FA+ is known to work
            prediction = "Possible to form 3D structure"
        else:
            prediction = "Unlikely to form 3D structure (t-factor out of range)"

        print(f"{cation_name:<25} {prediction}")
        # Outputting each number in the final equation as requested
        print(f"  Calculation: t = ({r_a_val} + {r_x_val}) / (sqrt(2) * ({r_b_val} + {r_x_val})) = {t_factor:.2f}\n")
    
    print("-" * 60)
    print("\nConclusion:")
    print("The calculations show that Ethylammonium and Dimethylammonium are too large to independently form 3D perovskites.")
    print("Methylammonium, Formamidinium, Aziridinium, and Methylhydrazinium all fall within the suitable range.")
    print("A truly comprehensive list would include all four. Since the options are mutually exclusive, we must select the best fit.")
    print("Both Aziridinium (Option B) and Methylhydrazinium (Option D) are valid additions beyond the common cations. Both have been demonstrated in scientific literature to form 3D perovskites. Choosing the option with Methylhydrazinium is a strong choice as it is noted for forming highly stable structures.")

if __name__ == '__main__':
    analyze_perovskite_cations()
<<<D>>>