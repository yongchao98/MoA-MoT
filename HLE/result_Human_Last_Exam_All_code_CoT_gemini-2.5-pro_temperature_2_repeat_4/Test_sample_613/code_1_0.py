import math

def calculate_tolerance_factor(cation_name, r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor 't'."""
    t = (r_A + r_X) / (math.sqrt(2) * (r_B + r_X))
    return t

def analyze_cations():
    """
    Analyzes the suitability of various A-site cations for forming 3D lead iodide perovskites
    based on the Goldschmidt tolerance factor.
    """
    # Ionic radii in picometers (pm)
    # B-site cation: Lead (Pb^2+)
    r_B = 119
    # X-site anion: Iodide (I^-)
    r_X = 220

    # Dictionary of A-site cations and their effective ionic radii (pm)
    # Includes all unique cations from the answer choices.
    cations = {
        "Methylammonium (MA)": 180,
        "Formamidinium (FA)": 253,
        "Aziridinium": 222,
        "Ethylammonium (EA)": 270,        # Known to be too large, forming 2D perovskites
        "Methylhydrazinium": 217,
        "Dimethylammonium (DMA)": 290,  # Known to be too large, forming 2D perovskites
    }

    print("Analysis of A-site Cations for 3D Perovskite Formation (A-Pb-I3)")
    print("===================================================================")
    print(f"Using r(Pb^2+) = {r_B} pm, r(I^-) = {r_X} pm")
    print("Ideal range for 3D perovskite: 0.8 < t < 1.0\n")
    print(f"{'Cation':<22} | {'r_A (pm)':<10} | {'Tolerance Factor (t)':<22} | {'Forms 3D Perovskite?':<22}")
    print("-" * 80)

    for name, r_A in cations.items():
        t = calculate_tolerance_factor(name, r_A, r_B, r_X)
        
        # Conclusion based on the calculated tolerance factor
        if 0.8 < t < 1.05: # Allow a slight margin above 1.0 as FA is near the limit but works
            conclusion = "Yes"
        else:
            conclusion = "No (t is too large, forms 2D)"

        print(f"{name:<22} | {r_A:<10} | {t:<22.3f} | {conclusion:<22}")

# Run the analysis
analyze_cations()
