def solve_structure():
    """
    This function analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    # 1. Molecular Formula Determination from Mass Spectrum
    mw = 135
    print(f"Step 1: Analysis of Mass Spectrum")
    print(f"The molecular ion peak (M+) is at m/z = {mw}.")
    print("An odd molecular weight suggests an odd number of Nitrogen atoms (Nitrogen Rule).")
    print("The molecular formula consistent with all data is C9H13N.")
    print(f"Mass of C9H13N = (9 * 12) + (13 * 1) + (1 * 14) = {9*12 + 13*1 + 14}, which matches the M+ peak.")
    print("Degrees of Unsaturation = 9 - (13/2) + (1/2) + 1 = 4. This corresponds to a benzene ring.\n")

    # 2. Fragment Identification from NMR data
    print("Step 2: Analysis of NMR Spectra")
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    print(f"13C NMR shows {len(c13_shifts)} signals: {c13_shifts}")
    print("DEPT-135 shows 1 negative signal (CH2), 5 positive signals (CH/CH3), and implies 1 quaternary C.")
    print("1H NMR shows signals for a phenyl group (5H), a CH group (1H), a CH2 group (2H), and a CH3 group (3H).")
    print("HSQC confirms the following C-H connections:")
    print("  - Proton at ~1.2 ppm (CH3) connects to Carbon at 19.2 ppm.")
    print("  - Protons at ~2.85 ppm (CH2) connect to Carbon at 43.5 ppm.")
    print("  - Proton at ~2.75 ppm (CH) connects to Carbon at 49.6 ppm.\n")

    # 3. Structure Assembly and Naming
    print("Step 3: Structure Elucidation and IUPAC Naming")
    print("The fragments identified are: a phenyl group (C6H5), a methylene group (-CH2-), a methine group (-CH-), a methyl group (-CH3), and an amine group (-NH2).")
    print("Assembling these fragments based on NMR data gives the structure: C6H5-CH2-CH(NH2)-CH3.\n")
    
    print("To determine the IUPAC name:")
    parent_chain_length = 3
    phenyl_position = 1
    amine_position = 2
    
    print(f"The longest carbon chain containing the principal functional group (amine) is {parent_chain_length} carbons long: 'propane'.")
    print(f"The amine group (-NH2) is at position number {amine_position}.")
    print(f"The phenyl group (-C6H5) is at position number {phenyl_position}.")
    print("\nCombining these parts, the final IUPAC name is:")
    print(f"({phenyl_position})-phenylpropan-({amine_position})-amine")


solve_structure()