def solve_structure():
    """
    Analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    print("Step-by-Step Analysis of the Spectroscopic Data:")
    print("================================================")

    # Step 1: Mass Spectrometry
    print("\n1. Mass Spectrum Analysis:")
    molecular_weight = 135
    print(f"The molecular ion (M+) peak is at m/z {molecular_weight}.")
    print(f"An odd molecular weight suggests an odd number of Nitrogen atoms. Assuming one N (mass=14).")
    print(f"Proposed molecular formula is C9H13N (Mass = 9*12 + 13*1 + 14 = 135).")

    # Step 2: Degree of Unsaturation
    print("\n2. Degree of Unsaturation (DoU):")
    C, H, N = 9, 13, 1
    dou = C - H/2 + N/2 + 1
    print(f"For C{C}H{H}N{N}, DoU = {C} - ({H}/2) + ({N}/2) + 1 = {int(dou)}.")
    print("A DoU of 4 indicates a benzene ring.")

    # Step 3: 13C and DEPT-135 NMR
    print("\n3. 13C NMR and DEPT-135 Analysis:")
    c13_shifts = {'ipso-C': 145.1, 'aromatic CH 1': 128.5, 'aromatic CH 2': 127.3, 'aromatic CH 3': 126.3,
                  'aliphatic CH': 49.6, 'aliphatic CH2': 43.5, 'aliphatic CH3': 19.2}
    print(f"The 13C shifts are: {list(c13_shifts.values())}.")
    print("DEPT-135 shows one negative signal (1 x CH2) and five positive signals (3 x aromatic CH, 1 x aliphatic CH, 1 x aliphatic CH3).")
    print("This confirms the side chain contains one CH2, one CH, and one CH3 group.")

    # Step 4: 1H NMR and HSQC
    print("\n4. 1H NMR and HSQC Analysis:")
    print("1H NMR shows a 5H multiplet (~7.2 ppm) for a C6H5- group.")
    print("It also shows a 3H doublet at ~1.2 ppm (CH3 next to a CH).")
    print("HSQC correlates the proton signals to their attached carbons:")
    print(f"  - CH proton at ~2.95 ppm correlates to Carbon at {c13_shifts['aliphatic CH']} ppm.")
    print(f"  - CH2 protons at ~2.85 ppm correlate to Carbon at {c13_shifts['aliphatic CH2']} ppm.")
    print(f"  - CH3 protons at ~1.2 ppm correlate to Carbon at {c13_shifts['aliphatic CH3']} ppm.")

    # Step 5: Structure Determination
    print("\n5. Structure Determination:")
    print("The fragments are C6H5-, -NH2, -CH-, -CH2-, and -CH3.")
    print("The connectivity is determined to be C6H5-CH2-CH(NH2)-CH3 based on the following:")
    print("  - The CH3 is a doublet, so it's next to the CH.")
    print("  - The chemical shift of the carbon attached to nitrogen (C-N) is expected to be higher than the benzylic carbon (C-Ph).")
    print(f"  - The data shows δC(CH) = {c13_shifts['aliphatic CH']} and δC(CH2) = {c13_shifts['aliphatic CH2']}.")
    print("  - This matches the structure where the CH is bonded to N and the CH2 is benzylic.")

    # Step 6: Final IUPAC Name
    iupac_name = "1-phenylpropan-2-amine"
    print("\n6. Conclusion:")
    print(f"The structure is confirmed as 1-phenylpropan-2-amine.")
    print(f"\nThe IUPAC name of the compound is: {iupac_name}")

solve_structure()