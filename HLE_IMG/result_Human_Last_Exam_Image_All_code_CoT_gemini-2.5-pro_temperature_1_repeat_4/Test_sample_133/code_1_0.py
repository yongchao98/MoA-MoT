def solve_mass_spec():
    """
    This function demonstrates the step-by-step deduction to identify the compound
    from its mass spectrum data.
    """
    # 1. Isotopic masses (nominal)
    mass_C = 12
    mass_H = 1
    mass_Cl35 = 35

    # 2. Key peaks from the spectrum data
    mz_molecular_ion_cluster_start = 260
    mz_base_peak = 225
    mz_fragment_2 = 190
    mz_fragment_3 = 155

    print("Step 1: Analyze the fragmentation pattern.")
    # 3. Calculate mass differences between major fragments
    loss_M_to_F1 = mz_molecular_ion_cluster_start - mz_base_peak
    loss_F1_to_F2 = mz_base_peak - mz_fragment_2
    loss_F2_to_F3 = mz_fragment_2 - mz_fragment_3

    print(f"The highest significant mass cluster starts at m/z = {mz_molecular_ion_cluster_start}.")
    print(f"The base peak (most intense fragment) is at m/z = {mz_base_peak}.")
    print(f"Mass lost from molecule to base peak: {mz_molecular_ion_cluster_start} - {mz_base_peak} = {loss_M_to_F1} amu.")
    print(f"This corresponds to the loss of one Chlorine atom (mass ~{mass_Cl35}).")
    
    print("\nThe fragmentation follows a clear cascade:")
    print(f"  {mz_base_peak} -> {mz_fragment_2} (Loss of {loss_F1_to_F2} amu, another Cl atom)")
    print(f"  {mz_fragment_2} -> {mz_fragment_3} (Loss of {loss_F2_to_F3} amu, a third Cl atom)")
    print("This indicates the molecule contains at least 3 Chlorine atoms and loses them sequentially.")

    print("\nStep 2: Determine the molecular formula.")
    print("The isotope pattern of the base peak at m/z 225 (with a large m/z 227 peak) indicates it contains 2 Chlorine atoms.")
    print("Since this fragment was formed by losing one Cl atom, the original molecule must have contained 3 Chlorine atoms.")

    # 4. Calculate the mass of the non-halogen part of the molecule ('R')
    mass_of_two_Cl = 2 * mass_Cl35
    mass_of_R = mz_base_peak - mass_of_two_Cl
    
    print(f"\nThe mass of the non-chlorine part (R) can be calculated from the base peak:")
    print(f"  mass(R) = m/z({mz_base_peak}) - mass(2 * Cl) = {mz_base_peak} - {mass_of_two_Cl} = {mass_of_R} amu.")

    # 5. Propose a formula for the 'R' group
    num_carbons_in_R = 12
    num_hydrogens_in_R = 11
    calculated_mass_R = num_carbons_in_R * mass_C + num_hydrogens_in_R * mass_H

    print(f"A plausible formula for a radical R with mass {mass_of_R} is C₁₂H₁₁.")
    print(f"  Calculation: {num_carbons_in_R} * {mass_C} (C) + {num_hydrogens_in_R} * {mass_H} (H) = {calculated_mass_R}")
    print(f"So the molecular formula is C₁₂H₁₁Cl₃.")

    print("\nStep 3: Propose the compound name.")
    print("Based on the formula C₁₂H₁₁Cl₃, a plausible structure with aromaticity is a derivative of phenylcyclohexene.")
    final_answer = "3,4,5-trichloro-1-phenylcyclohex-1-ene"
    print(f"\nFinal identified compound (IUPAC Name): {final_answer}")

solve_mass_spec()