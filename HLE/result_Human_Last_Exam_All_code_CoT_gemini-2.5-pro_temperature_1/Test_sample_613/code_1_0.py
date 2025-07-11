import math

def calculate_and_explain_tolerance_factor():
    """
    Calculates and explains the Goldschmidt tolerance factor for various A-site cations
    in a lead bromide perovskite structure (A-Pb-Br3).
    """

    # Define ionic radii in picometers (pm). Effective radii are used for molecules.
    radii = {
        "A_cations": {
            "Cesium (Cs+)": 167,
            "Methylammonium (MA+)": 217,
            "Formamidinium (FA+)": 253,
            "Aziridinium (AZ+)": 202,
            "Ethylammonium (EA+)": 274,
            "Methylhydrazinium (MHy+)": 217,
            "Dimethylammonium (DMA+)": 272
        },
        "B_cation": {
            "Lead (Pb2+)": 119
        },
        "X_anion": {
            "Bromide (Br-)": 196
        }
    }

    r_B = radii["B_cation"]["Lead (Pb2+)"]
    r_X = radii["X_anion"]["Bromide (Br-)"]
    sqrt2 = math.sqrt(2)

    print("Goldschmidt Tolerance Factor (t) Analysis for A-Pb-Br3 Perovskites")
    print("Formula: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
    print("A stable 3D perovskite structure is favored when 't' is in the range of 0.8 to 1.0.")
    print("-" * 75)

    for name, r_A in radii["A_cations"].items():
        t = (r_A + r_X) / (sqrt2 * (r_B + r_X))
        
        # Determine suitability
        if t > 1.05: # A little tolerance above 1.0, but >1.05 is almost certainly not 3D
            suitability = "No (Too large, typically forms 2D structures)"
        elif t < 0.8:
            suitability = "No (Too small, favors non-perovskite phases)"
        else:
            suitability = "Yes (Geometrically suitable for 3D structure)"

        print(f"For {name}:")
        # Print the full equation with numbers as requested
        equation_str = f"t = ({r_A} + {r_X}) / (1.414 * ({r_B} + {r_X})) = {r_A + r_X} / (1.414 * {r_B + r_X}) = {t:.3f}"
        print(equation_str)
        print(f"Result: {suitability}\n")

    print("--- Conclusion ---")
    print("The calculations show that Ethylammonium (t=1.055) and Dimethylammonium (t=1.050) are too large to fit in the 3D perovskite structure, which rules out answer choices C and E.")
    print("Cesium, Methylammonium, Formamidinium, Aziridinium, and Methylhydrazinium all have tolerance factors within the acceptable range.")
    print("Therefore, the most comprehensive list among the choices that includes only cations capable of forming 3D structures is the one containing Cesium, Methylammonium, Formamidinium, and Methylhydrazinium.")

calculate_and_explain_tolerance_factor()