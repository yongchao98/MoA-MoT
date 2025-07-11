def solve_structure():
    """
    This function provides a step-by-step analysis of the provided spectroscopic data
    to determine the IUPAC name of the unknown compound.
    """
    print("--- Spectroscopic Analysis ---")
    
    # Step 1: Mass Spectrometry
    mw = 135
    formula = "C9H13N"
    dbe = 4
    print(f"\nStep 1: Mass Spectrum Analysis")
    print(f"The molecular ion peak (M+) at m/z {mw} suggests an odd molecular weight.")
    print(f"By the Nitrogen Rule, this indicates an odd number of nitrogen atoms (likely one).")
    print(f"The molecular formula is determined to be {formula}.")
    print(f"The Degree of Unsaturation is {dbe}, consistent with a benzene ring.")

    # Step 2: 13C NMR and DEPT-135
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    print(f"\nStep 2: 13C NMR and DEPT-135 Analysis")
    print(f"13C NMR signals (ppm): {c13_shifts}")
    print("DEPT-135 data indicates 1 negative signal (CH2) and 5 positive signals (CH/CH3).")
    print("Combined with the 13C shifts, we can deduce the following carbon types:")
    print(" - Aromatic: 1 quaternary C (~145.1), 3 CH groups (~126-129)")
    print(" - Aliphatic: 1 CH2 group, 1 CH group, 1 CH3 group")

    # Step 3: 1H NMR and HSQC
    print("\nStep 3: 1H NMR and HSQC Analysis")
    print("1H NMR shows signals for:")
    print(" - A 5H multiplet (~7.2 ppm): Monosubstituted phenyl ring (C6H5-)")
    print(" - A 2H multiplet (~2.7 ppm): A CH2 group")
    print(" - A 1H multiplet (~2.9 ppm): A CH group")
    print(" - A 3H doublet (~1.2 ppm): A CH3 group next to a CH")
    print("HSQC data correlates these protons to their respective carbons:")
    print(" - CH3 (1H ~1.2 ppm) <=> C (~19.2 ppm)")
    print(" - CH2 (1H ~2.7 ppm) <=> C (~49.6 ppm)")
    print(" - CH (1H ~2.9 ppm) <=> C (~43.5 ppm)")

    # Step 4: Structure Assembly and Verification
    print("\nStep 4: Structure Assembly and Verification")
    print("The fragments are: C6H5-, -CH-, -CH2-, -CH3, and -NH2.")
    print("Assembling these fragments leads to the structure: C6H5-CH(CH3)-CH2-NH2.")
    print("This structure is verified by the MS fragmentation pattern:")
    print(" - The base peak at m/z 30 corresponds to the [CH2=NH2]+ fragment, formed by alpha-cleavage.")
    print(" - The peak at m/z 91 corresponds to the [C7H7]+ (tropylium) ion, characteristic of a benzyl group.")

    # Step 5: IUPAC Name
    iupac_name = "1-phenylpropan-2-amine"
    print("\nStep 5: Final IUPAC Name")
    print(f"The parent chain is propane. An amino group is on carbon 2, and a phenyl group is on carbon 1.")
    print(f"The final IUPAC name is: {iupac_name}")

solve_structure()