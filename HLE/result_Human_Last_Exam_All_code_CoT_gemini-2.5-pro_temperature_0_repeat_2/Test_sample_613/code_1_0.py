import math

def calculate_tolerance_factor():
    """
    Calculates and evaluates the Goldschmidt tolerance factor for various A-site
    cations in a lead bromide perovskite structure (APbBr3).
    """
    # Ionic and effective radii in picometers (pm)
    radii = {
        'Pb2+': 119,
        'Br-': 196,
        'Cs+': 167,
        'Methylammonium (MA+)': 180,
        'Formamidinium (FA+)': 253,
        'Ethylammonium (EA+)': 230,
        'Dimethylammonium (DMA+)': 272,
        'Methylhydrazinium (MHy+)': 217,
        'Aziridinium (AZ+)': 202,
    }

    r_B = radii['Pb2+']
    r_X = radii['Br-']
    sqrt_2 = math.sqrt(2)

    print("Evaluating A-site cations for APbBr3 perovskite formation.")
    print(f"Using Goldschmidt Tolerance Factor: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print(f"For APbBr3: r_B(Pb2+) = {r_B} pm, r_X(Br-) = {r_X} pm")
    print("A stable 3D structure generally requires 0.8 < t < 1.0.\n")

    # Cations to evaluate from the answer choices
    cations_to_check = [
        'Cs+',
        'Methylammonium (MA+)',
        'Formamidinium (FA+)',
        'Ethylammonium (EA+)',
        'Dimethylammonium (DMA+)',
        'Methylhydrazinium (MHy+)',
        'Aziridinium (AZ+)', # Included for completeness
    ]

    for cation_name in sorted(cations_to_check):
        r_A = radii[cation_name]
        
        numerator = r_A + r_X
        denominator = sqrt_2 * (r_B + r_X)
        t = numerator / denominator
        
        print(f"--- Analysis for {cation_name} ---")
        # Print the full equation with substituted numbers
        print(f"Equation: t = ({r_A} + {r_X}) / (sqrt(2) * ({r_B} + {r_X}))")
        print(f"Calculation: t = {numerator} / {denominator:.1f} = {t:.3f}")

        # Analysis based on 't' value and experimental facts
        if cation_name in ['Cs+', 'Methylammonium (MA+)', 'Formamidinium (FA+)']:
            print("Conclusion: Forms stable 3D perovskites. A core A-site cation.")
        elif cation_name == 'Methylhydrazinium (MHy+)':
            print("Conclusion: 't' is in the ideal range, and it is experimentally confirmed to form a stable 3D perovskite with bromide.")
        elif cation_name == 'Ethylammonium (EA+)':
            print("Conclusion: 't' is in the acceptable range, but the cation's shape promotes the formation of 2D layered structures, not 3D.")
        elif cation_name == 'Dimethylammonium (DMA+)':
            print("Conclusion: 't' is too large (> 1.0). This cation is too bulky to form a 3D perovskite structure.")
        elif cation_name == 'Aziridinium (AZ+)':
            print("Conclusion: 't' is in the ideal range, but experimentally it forms a 1D structure with bromide, though it can form a 3D structure with the larger iodide anion.")
        print("-" * 40 + "\n")

    print("Final Verdict:")
    print("The most comprehensive list of cations that can INDEPENDENTLY form 3D lead halide perovskites includes the standard Cesium, Methylammonium, and Formamidinium, plus Methylhydrazinium, which is experimentally proven to form a stable 3D structure with the example halide (bromide). Ethylammonium and Dimethylammonium are not suitable for forming independent 3D structures.")

if __name__ == '__main__':
    calculate_tolerance_factor()