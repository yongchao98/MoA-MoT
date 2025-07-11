def solve_structure():
    """
    Analyzes the provided spectral data step-by-step to determine the IUPAC name of the compound.
    """
    
    print("Step 1: Analysis of Mass Spectrum (MS)")
    molecular_ion = 135
    print(f"The mass spectrum shows a molecular ion peak (M+) at m/z = {molecular_ion}.")
    print("According to the Nitrogen Rule, an odd molecular weight suggests the presence of an odd number of nitrogen atoms.")
    print("Let's assume one nitrogen atom (N, atomic weight ~14).")
    remaining_mass = molecular_ion - 14
    print(f"Remaining mass for Carbon and Hydrogen = {molecular_ion} - 14 = {remaining_mass}.")
    num_carbons = 9
    num_hydrogens = 13
    print(f"A plausible molecular formula is C{num_carbons}H{num_hydrogens}N, which has a mass of (9*12 + 13*1 + 14) = {num_carbons*12 + num_hydrogens*1 + 14}.")
    # DoU = C + 1 - H/2 + N/2
    degree_of_unsaturation = num_carbons + 1 - (num_hydrogens / 2) + (1 / 2)
    print(f"The degree of unsaturation for C9H13N is {int(degree_of_unsaturation)}, suggesting an aromatic ring (DoU=4).\n")
    
    print("Step 2: Analysis of Infrared (IR) Spectrum")
    print("- Peaks just above 3000 cm-1: Aromatic C-H stretch.")
    print("- Peaks just below 3000 cm-1: Aliphatic C-H stretch.")
    print("- Sharp peak around 3400 cm-1: N-H stretch (likely from a primary or secondary amine).")
    print("- Peaks around 1600 & 1500 cm-1: C=C stretch, characteristic of an aromatic ring.")
    print("Conclusion: The molecule contains an aromatic ring, aliphatic groups, and an N-H bond.\n")
    
    print("Step 3: Analysis of 13C NMR and DEPT-135 Spectra")
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    print(f"There are {len(c13_shifts)} unique carbon signals, which means there is symmetry in the C9 molecule.")
    print("- Four signals are in the aromatic region (120-150 ppm), consistent with a monosubstituted benzene ring (C-ipso, C-ortho, C-meta, C-para).")
    print("- Three signals are in the aliphatic region (< 60 ppm): 49.6, 43.5, 19.2.")
    print("The DEPT-135 spectrum shows one negative signal (CH2) and five positive signals (CH and CH3).")
    print("- The aromatic CH groups account for 3 positive signals.")
    print("- Therefore, the 3 aliphatic carbons must consist of one CH2 group (negative), one CH group (positive), and one CH3 group (positive).")
    print("Carbon Skeleton: A monosubstituted phenyl group (C6H5-) and a C3 side chain containing one CH2, one CH, and one CH3 group.\n")

    print("Step 4: Analysis of 1H NMR and HSQC Spectra")
    print("- Aromatic region (~7.2 ppm): A multiplet integrating to 5H, confirming the C6H5- (phenyl) group.")
    print("- Aliphatic region:")
    print("  - A multiplet at ~2.8-2.9 ppm (2H), correlated by HSQC to the carbon at 43.5 ppm. This is the CH2 group.")
    print("  - A multiplet at ~2.6-2.7 ppm (1H), correlated by HSQC to the carbon at 49.6 ppm. This is the CH group.")
    print("  - A doublet at ~1.2 ppm (3H), correlated by HSQC to the carbon at 19.2 ppm. This is the CH3 group.")
    print("The total proton count is 5 (Aromatic) + 2 (CH2) + 1 (CH) + 3 (CH3) + 2 (from NH2) = 13H. This matches the formula C9H13N.\n")

    print("Step 5: Proposing a Structure")
    print("Based on the fragments (C6H5-, -CH2-, -CH-, -CH3, and an -NH2 group) and their connectivity from NMR splitting patterns:")
    print("- The 3H doublet at 1.2 ppm (-CH3) is coupled to a single proton, which must be the -CH- group.")
    print("- The 1H multiplet (-CH-) must be coupled to the -CH3 (3H) and the -CH2- (2H), for a total of 5 neighbors (expecting a sextet, observed as a multiplet).")
    print("- The 2H multiplet (-CH2-) is coupled to the -CH- group (1H), which would suggest a doublet, but appears more complex, possibly due to being adjacent to a stereocenter.")
    print("Assembling these pieces leads to the only possible structure: Phenyl-CH2-CH(NH2)-CH3.\n")
    
    print("Step 6: Final Structure Verification")
    print("Proposed Structure: 1-phenylpropan-2-amine")
    print("- Formula: C9H13N. Correct.")
    print("- Molecular Weight: 135. Correct.")
    print("- IR: Has aromatic, aliphatic, and N-H bonds. Correct.")
    print("- 13C NMR: Predicts 7 signals (4 aromatic, 3 aliphatic). Correct.")
    print("- DEPT-135: Predicts 1 CH2 (negative) and 5 positive (4 CH, 1 CH3). Correct.")
    print("- 1H NMR: Predicts C6H5 (5H), CH2 (2H), CH (1H), CH3 (3H) with the observed splitting. Correct.")
    print("- MS fragmentation: A peak at m/z 91 (tropylium ion from benzyl cleavage) is expected and observed.\n")

    print("Step 7: Final Answer - IUPAC Name")
    final_name = "1-phenylpropan-2-amine"
    print(f"The IUPAC name of the compound is: {final_name}")

solve_structure()
<<<1-phenylpropan-2-amine>>>