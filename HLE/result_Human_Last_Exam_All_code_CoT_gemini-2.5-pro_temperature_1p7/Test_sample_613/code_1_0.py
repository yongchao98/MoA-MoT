import math

def calculate_goldschmidt_tolerance_factor():
    """
    Calculates the Goldschmidt tolerance factor for various A-site cations in a 
    lead bromide (A-Pb-Br3) perovskite structure to determine their suitability
    for forming a stable 3D lattice.
    """
    # Ionic radii in picometers (pm). Sources for organic cations can vary slightly.
    # These are commonly accepted effective radii.
    radii = {
        'A_cations': {
            'Methylammonium (MA)': 217,
            'Formamidinium (FA)': 253,
            'Ethylammonium (EA)': 274,
            'Dimethylammonium (DMA)': 272,
            'Methylhydrazinium (MHy)': 217,
            'Aziridinium (AZ)': 258,
        },
        'B_cation': {
            'Lead (Pb2+)': 119
        },
        'X_anion': {
            'Bromide (Br-)': 196
        }
    }

    r_B = radii['B_cation']['Lead (Pb2+)']
    r_X = radii['X_anion']['Bromide (Br-)']
    
    print("Calculating Goldschmidt Tolerance Factor (t) for A-Pb-Br3 perovskites.")
    print("Formula: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print(f"Using radii: r(Pb2+) = {r_B} pm, r(Br-) = {r_X} pm\n")
    print("-" * 60)

    denominator = math.sqrt(2) * (r_B + r_X)

    # List of all unique cations from the options
    all_cations = [
        'Methylammonium (MA)', 
        'Formamidinium (FA)',
        'Aziridinium (AZ)', 
        'Ethylammonium (EA)', 
        'Methylhydrazinium (MHy)',
        'Dimethylammonium (DMA)'
    ]
    
    # Analyze each cation from the answer choices
    for cation_name in sorted(all_cations, key=lambda x: radii['A_cations'][x]):
        r_A = radii['A_cations'][cation_name]
        
        numerator = r_A + r_X
        t = numerator / denominator

        # Determine suitability based on t
        if 0.8 < t < 1.0:
            conclusion = "Forms stable 3D perovskite."
        elif t >= 1.0:
            conclusion = "Too large. Tends to form 2D or distorted structures."
        else: # t <= 0.8
            conclusion = "Too small. Tends to form non-perovskite structures."
            
        print(f"Cation: {cation_name}")
        print(f"Equation: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X}))")
        print(f"Calculation: t = {numerator} / {denominator:.1f} = {t:.3f}")
        print(f"Conclusion: {conclusion}")
        print("-" * 60)

# Execute the calculation
calculate_goldschmidt_tolerance_factor()