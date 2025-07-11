def identify_compound():
    """
    Analyzes spectroscopic data to determine the IUPAC name of a compound.
    """
    # --- Data from the problem ---
    molecular_ion_mz = 135
    c13_shifts_ppm = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    h1_shifts_ppm = {"C6H5": "7.2-7.4", "CH2": "~2.8", "CH": "~2.7", "CH3": "~1.2", "NH2": "~1.0"}
    hsqc_correlations = {
        "CH3": ("~1.2", 19.2),
        "CH": ("~2.7", 49.6),
        "CH2": ("~2.8", 43.5)
    }

    # --- Analysis ---
    print("Step-by-step analysis of spectroscopic data:")
    
    # 1. Mass Spectrometry
    print("\n1. Mass Spectrum (MS) Analysis:")
    print(f"The molecular ion (M+) peak is at m/z = {molecular_ion_mz}.")
    print("The odd molecular weight implies the presence of one nitrogen atom (Nitrogen Rule).")
    print("The molecular formula is C9H13N.")
    # DBE = C + 1 - H/2 + N/2 = 9 + 1 - 13/2 + 1/2 = 4
    dbe = 4
    print(f"The Degree of Unsaturation is {dbe}, suggesting a benzene ring.")

    # 2. IR and NMR Spectroscopy
    print("\n2. IR, 13C NMR, DEPT-135, and 1H NMR Analysis:")
    print("- IR spectrum shows a primary amine (-NH2), a benzene ring, and aliphatic C-H bonds.")
    print(f"- 13C NMR shows 7 signals: {c13_shifts_ppm}, consistent with a monosubstituted benzene ring and a 3-carbon side chain.")
    print("- DEPT-135 indicates the side chain contains one CH2, one CH, and one CH3 group.")
    print(f"- 1H NMR confirms the C6H5 group ({h1_shifts_ppm['C6H5']} ppm), a -CH-CH3 fragment (doublet at {h1_shifts_ppm['CH3']} ppm), a CH2 group ({h1_shifts_ppm['CH2']} ppm), a CH group ({h1_shifts_ppm['CH']} ppm), and the NH2 group ({h1_shifts_ppm['NH2']} ppm).")

    # 3. HSQC and Structure Determination
    print("\n3. HSQC Analysis and Final Structure Determination:")
    print("HSQC confirms the following direct C-H bonds:")
    print(f"  - CH3 group: 1H at {hsqc_correlations['CH3'][0]} ppm is attached to 13C at {hsqc_correlations['CH3'][1]} ppm.")
    print(f"  - CH group: 1H at {hsqc_correlations['CH'][0]} ppm is attached to 13C at {hsqc_correlations['CH'][1]} ppm.")
    print(f"  - CH2 group: 1H at {hsqc_correlations['CH2'][0]} ppm is attached to 13C at {hsqc_correlations['CH2'][1]} ppm.")
    print("\nComparing these values to expected shifts for possible isomers:")
    print("The structure is identified as 1-phenylpropan-2-amine (C6H5-CH2-CH(NH2)-CH3).")
    print(f"  - The benzylic CH2 group's 13C shift ({hsqc_correlations['CH2'][1]} ppm) matches the expected value (~40-45 ppm).")
    print(f"  - The CH group attached to the amine's 13C shift ({hsqc_correlations['CH'][1]} ppm) matches the expected value (~50-60 ppm).")
    print("This assignment is consistent across all provided spectral data.")

    # 4. Final Answer
    print("\n-------------------------------------------------")
    print("The IUPAC name of the compound is:")
    print("1-phenylpropan-2-amine")
    print("-------------------------------------------------")

if __name__ == '__main__':
    identify_compound()