import collections

def solve_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the unknown compound.
    """
    # Mass spectrum data from the table (key fragments)
    spectrum = {
        266: 3.22, 264: 12.9, 262: 29.0, 260: 34.8,
        227: 66.7, 225: 100.0,
        192: 24.5, 190: 48.4,
        155: 16.1,
        143: 27.8
    }
    
    # Atomic masses (integer mass for nominal calculation)
    atomic_mass = {'H': 1, 'C': 12, 'N': 14, 'O': 16, 'Cl': 35, '37Cl': 37}

    print("Step 1 & 2: Analyze Molecular Ion and Nitrogen Rule")
    # The highest m/z cluster is at 260, 262, 264, 266.
    # The peak for the molecule with the lightest isotopes (M) is at m/z 260.
    M_peak_mz = 260
    print(f"The molecular ion peak (M+) with the lightest isotopes appears to be at m/z = {M_peak_mz}.")
    # Nitrogen Rule: An even molecular weight implies an even number of nitrogen atoms (usually 0 or 2).
    print(f"The nominal mass ({M_peak_mz}) is even, which suggests the compound has 0 or an even number of nitrogen atoms.\n")

    print("Step 3: Analyze Isotope Patterns")
    # The M+ cluster has peaks at M, M+2, M+4, M+6. This pattern is characteristic of a compound with three chlorine atoms.
    # Theoretical ratio for 3 Cl atoms is approx. 100 : 96 : 31 : 3.
    # Observed ratio (normalized to m/z 260):
    I_260 = spectrum[260]
    ratio_262 = (spectrum[262] / I_260) * 100
    ratio_264 = (spectrum[264] / I_260) * 100
    ratio_266 = (spectrum[266] / I_260) * 100
    print(f"The M+ cluster at m/z {M_peak_mz}, {M_peak_mz+2}, {M_peak_mz+4}, {M_peak_mz+6} suggests 3 chlorine atoms.")
    print(f"Observed intensity ratio (normalized): 100 : {ratio_262:.1f} : {ratio_264:.1f} : {ratio_266:.1f}")
    print("This pattern confirms the presence of 3 Cl atoms, despite some deviation from theoretical values.\n")
    
    print("Step 4 & 5: Analyze Base Peak and Fragmentation")
    base_peak_mz = 225
    print(f"The base peak (most intense fragment) is at m/z = {base_peak_mz}.")
    
    # Check the isotope pattern of the base peak
    base_peak_Mplus2_mz = 227
    base_peak_ratio = (spectrum[base_peak_Mplus2_mz] / spectrum[base_peak_mz]) * 100
    print(f"The base peak has an isotope peak at m/z {base_peak_Mplus2_mz} with relative intensity {base_peak_ratio:.1f}%.")
    print("A 100:67 ratio is characteristic of a fragment containing 2 chlorine atoms.\n")

    print("Analyzing the fragmentation pathway:")
    loss_to_base_peak = M_peak_mz - base_peak_mz
    print(f"Fragment 1: Loss from M+ ({M_peak_mz}) to base peak ({base_peak_mz}) is {loss_to_base_peak} amu.")
    print(f"This corresponds to the loss of one chlorine atom ({atomic_mass['Cl']} amu).")
    print(f"So, the base peak [M-Cl]+ has a formula of [C12H11Cl2]+.\n")

    frag2_mz = 190
    loss_from_base_peak = base_peak_mz - frag2_mz
    print(f"Fragment 2: Loss from base peak ({base_peak_mz}) to fragment at m/z {frag2_mz} is {loss_from_base_peak} amu.")
    print("This corresponds to the loss of a second chlorine atom.")
    print(f"This fragment [M-2Cl]+ has a formula of [C12H11Cl]+.\n")
    
    frag3_mz = 155
    loss_from_frag2 = frag2_mz - frag3_mz
    print(f"Fragment 3: Loss from fragment at m/z {frag2_mz} to fragment at m/z {frag3_mz} is {loss_from_frag2} amu.")
    print("This corresponds to the loss of a third chlorine atom.")
    print(f"This fragment [M-3Cl]+ has a formula of [C12H11]+.\n")
    print("This sequential loss of three Cl atoms strongly suggests a -CCl3 (trichloromethyl) group.\n")

    print("Step 6: Deduce Molecular Formula")
    # From M=260 and 3 Cl atoms, we can find the mass of the rest of the molecule.
    mass_of_Cl3 = 3 * atomic_mass['Cl']
    mass_of_hydrocarbon = M_peak_mz - mass_of_Cl3
    print(f"Mass of three Cl atoms = 3 * {atomic_mass['Cl']} = {mass_of_Cl3}")
    print(f"Mass of the rest of the molecule (CxHy) = {M_peak_mz} - {mass_of_Cl3} = {mass_of_hydrocarbon}")
    # Find formula for CxHy = 155
    num_C = mass_of_hydrocarbon // atomic_mass['C']
    num_H = mass_of_hydrocarbon % atomic_mass['C']
    print(f"A hydrocarbon part with mass 155 corresponds to C{num_C}H{num_H}.")
    print(f"Molecular Formula = C{num_C}H{num_H}Cl3.\n")
    
    print("Step 7: Propose Structure and Name")
    print("The fragmentation suggests a structure of R-CCl3, where R is C11H11.")
    # The cleavage of the C-CCl3 bond would produce an [R]+ fragment.
    mass_R = (num_C - 1) * atomic_mass['C'] + num_H * atomic_mass['H']
    print(f"This predicts a fragment [C{num_C-1}H{num_H}]+ at m/z = {mass_R}.")
    print(f"A significant peak is observed in the spectrum at m/z = {143}, confirming this structural feature.")
    print("\nA plausible structure for C12H11Cl3 with these features is 1-phenyl-2-(trichloromethyl)cyclopent-1-ene.")
    print("This structure is consistent with the formula and all major fragments observed.\n")

    print("Final Answer:")
    final_name = "1-phenyl-2-(trichloromethyl)cyclopent-1-ene"
    print(f"The compound is {final_name}.")

solve_mass_spectrum()