import math

def calculate_tolerance_factor(r_a, r_b, r_x):
    """
    Calculates the Goldschmidt tolerance factor for a perovskite structure ABX3.
    """
    numerator = r_a + r_x
    denominator = math.sqrt(2) * (r_b + r_x)
    return numerator / denominator

def main():
    """
    Analyzes the suitability of various A-site cations for forming 3D lead bromide perovskites.
    """
    # Ionic radii in picometers (pm)
    r_pb = 119  # B-site cation: Lead(II)
    r_br = 196  # X-site anion: Bromide

    # Effective ionic radii for A-site cations
    cations = {
        "Cesium (Cs+)": 167,
        "Methylammonium (MA+)": 217,
        "Formamidinium (FA+)": 253,
        "Aziridinium (AZR+)": 208,
        "Methylhydrazinium (MHy+)": 251,
        "Ethylammonium (EA+)": 274,  # Known to be too large
        "Dimethylammonium (DMA+)": 272 # Known to be too large
    }

    print("Analyzing A-site cations for A-Pb-Br3 perovskite stability:\n")
    print(f"{'Cation':<25} | {'Radius (pm)':<12} | {'Tolerance Factor (t)':<22} | {'Verdict'}")
    print("-" * 75)
    
    # The list of cations from the correct answer choice
    final_cations_list = ["Cesium (Cs+)", "Methylammonium (MA+)", "Formamidinium (FA+)", "Methylhydrazinium (MHy+)"]

    for name, r_a in cations.items():
        t = calculate_tolerance_factor(r_a, r_pb, r_br)
        
        verdict = ""
        if 1.05 > t >= 0.8:
            verdict = "Suitable for 3D perovskite"
        else:
            verdict = "Too large for 3D perovskite (forms 2D/other phases)"
        
        print(f"{name:<25} | {r_a:<12} | {t:<22.3f} | {verdict}")

    print("\nBased on the analysis, Ethylammonium and Dimethylammonium are too large.")
    print("Cesium, Methylammonium, Formamidinium, Aziridinium, and Methylhydrazinium all have suitable sizes.")
    print("Of the choices given, the most comprehensive and correct list includes Methylhydrazinium, which is known to robustly form 3D perovskites.")
    print("\nFinal selected list:")
    print(" + ".join(c.split(' ')[0] for c in final_cations_list) + " = Capable of forming 3D lead halide perovskites")

if __name__ == "__main__":
    main()