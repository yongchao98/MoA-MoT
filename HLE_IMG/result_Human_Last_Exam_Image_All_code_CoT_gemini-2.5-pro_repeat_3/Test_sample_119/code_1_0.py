def solve_structure():
    """
    This function analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    
    # --- Data Summary ---
    mw = 135
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    dept_135 = {"negative": 1, "positive": 5}
    h1_shifts = {
        "aromatic": "~7.2 ppm (5H, multiplet)",
        "aliphatic_ch2": "~2.8 ppm (2H, multiplet)",
        "aliphatic_ch": "~2.6 ppm (1H, multiplet)",
        "aliphatic_ch3": "~1.2 ppm (3H, doublet)",
        "amine_nh2": "~1.1 ppm (2H, broad)"
    }
    hsqc_correlations = {
        "1.2 ppm H": "19.2 ppm C",
        "2.6 ppm H": "49.6 ppm C",
        "2.8 ppm H": "43.5 ppm C",
        "7.2 ppm H": "126.3, 127.3, 128.5 ppm C"
    }
    
    # --- Analysis ---
    print("Step-by-step Analysis:")
    
    print("\n1. Molecular Formula and Degree of Unsaturation (DoU):")
    print(f"   - The mass spectrum shows a molecular ion peak at m/z = {mw}.")
    print("   - An odd molecular weight suggests one nitrogen atom.")
    print("   - The IR spectrum shows N-H stretches (primary amine) and aromatic C-H/C=C stretches.")
    print("   - The proposed molecular formula is C9H13N (Calculated MW = 9*12 + 13*1 + 14 = 135).")
    dou = 9 + 1 - (13 / 2) + (1 / 2)
    print(f"   - The Degree of Unsaturation is {int(dou)}, which corresponds to a benzene ring.")

    print("\n2. Carbon Skeleton from 13C and DEPT-135 NMR:")
    print(f"   - There are {len(c13_shifts)} unique carbon signals.")
    print(f"   - DEPT-135 shows {dept_135['negative']} negative signal (1 x CH2) and {dept_135['positive']} positive signals (CHs and CH3s).")
    print(f"   - The signal at {c13_shifts[0]} ppm is a quaternary carbon (no H attached).")
    print(f"   - The signals at {c13_shifts[1]}, {c13_shifts[2]}, {c13_shifts[3]} ppm are aromatic CH carbons.")
    print(f"   - The remaining signals {c13_shifts[4]}, {c13_shifts[5]}, {c13_shifts[6]} must be one CH, one CH2, and one CH3.")

    print("\n3. Fragment Assembly from 1H NMR and HSQC:")
    print(f"   - 1H NMR shows a monosubstituted phenyl group ({h1_shifts['aromatic']}).")
    print(f"   - A methyl group ({h1_shifts['aliphatic_ch3']}) is a doublet, so it's next to a CH group.")
    print(f"   - This gives the fragment: -CH(CH3). The remaining fragment is -CH2-.")
    print(f"   - The structure must be C6H5-CH2-CH(NH2)-CH3.")

    print("\n4. Final Verification with HSQC:")
    print(f"   - The CH3 protons ({list(hsqc_correlations.keys())[0]}) connect to the C at {hsqc_correlations['1.2 ppm H']}.")
    print(f"   - The CH proton ({list(hsqc_correlations.keys())[1]}) connects to the C at {hsqc_correlations['2.6 ppm H']}.")
    print(f"   - The CH2 protons ({list(hsqc_correlations.keys())[2]}) connect to the C at {hsqc_correlations['2.8 ppm H']}.")
    print("   - This confirms the structure.")

    print("\n--- Conclusion ---")
    final_name = "1-phenylpropan-2-amine"
    print(f"The IUPAC name of the compound is: {final_name}")

solve_structure()