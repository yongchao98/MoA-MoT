def identify_compound_and_name():
    """
    This script outlines the analysis of spectroscopic data to identify a chemical compound and provides its IUPAC name.
    """
    print("Analysis of Spectroscopic Data")
    print("==============================\n")

    # 1. Mass Spectrometry Analysis
    print("1. Mass Spectrum (MS):")
    molecular_weight = 135
    formula = "C9H13N"
    print(f"- Molecular ion peak (M+) at m/z = {molecular_weight}, suggesting the molecular weight is {molecular_weight} u.")
    print(f"- The odd molecular weight points to the presence of one nitrogen atom. The molecular formula is determined to be {formula}.")
    print("- Degrees of Unsaturation = 9 + 1 - (13/2) + (1/2) = 4, which indicates a benzene ring.\n")

    # 2. Infrared Spectroscopy Analysis
    print("2. Infrared (IR) Spectrum:")
    print("- A doublet near 3300-3400 cm^-1 indicates a primary amine (-NH2).")
    print("- Peaks for aromatic C-H (>3000 cm^-1) and aliphatic C-H (<3000 cm^-1) are present.")
    print("- Strong absorptions around 700-750 cm^-1 suggest a monosubstituted benzene ring.\n")

    # 3. NMR Spectroscopy Analysis
    print("3. 13C NMR, DEPT-135, and 1H NMR Analysis:")
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    print(f"- 13C NMR signals at {c13_shifts} ppm.")
    print("- DEPT-135 shows one negative signal (CH2) and five positive signals (CH and CH3).")
    print("- This means the aliphatic side chain contains one CH2, one CH, and one CH3 group.")
    print("- 1H NMR shows a 5H multiplet for a C6H5- group, a 3H doublet (~1.2 ppm), a 2H multiplet (~2.8 ppm), and a 1H multiplet (~2.6-3.0 ppm).")
    print("- The 3H doublet and 1H multiplet confirm a -CH-CH3 fragment.\n")

    # 4. Structure Elucidation
    print("4. Structure Determination:")
    print("Combining the fragments (C6H5-, -CH2-, -CH-, -CH3, -NH2) and NMR coupling information leads to the structure: C6H5-CH2-CH(NH2)-CH3.")
    print("HSQC data confirms the following connections:")
    print("  - H(~1.2 ppm) <=> C(19.2 ppm) : CH3 group")
    print("  - H(~2.8 ppm) <=> C(43.5 ppm) : CH2 group")
    print("  - H(~2.6-3.0 ppm) <=> C(49.6 ppm) : CH group\n")
    
    # 5. IUPAC Nomenclature
    print("5. IUPAC Name:")
    iupac_name = "1-phenylpropan-2-amine"
    print(f"The longest carbon chain with the amine is a 3-carbon chain (propane).")
    print(f"The amine is at position 2, and the phenyl group is at position 1.")
    print(f"The final IUPAC name is: {iupac_name}")

identify_compound_and_name()