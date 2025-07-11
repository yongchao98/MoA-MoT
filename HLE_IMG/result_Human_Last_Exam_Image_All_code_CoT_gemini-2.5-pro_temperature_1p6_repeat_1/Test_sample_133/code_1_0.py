def solve_mass_spectrometry():
    """
    Solves the mass spectrometry problem by calculating and explaining the key fragments.
    """
    # Key masses from the spectrum table
    mass_molecular_ion = 262  # M+ peak for the molecule with 2 * 35Cl
    mass_base_peak = 225      # The most intense peak
    mass_M_minus_Cl = 227     # Peak corresponding to the loss of one 35Cl
    
    # Atomic masses for calculation
    mass_Cl35 = 35
    mass_H = 1
    mass_C = 12

    # Step 1: Analyze the molecular ion and determine the non-halogen mass
    mass_of_two_Cl = 2 * mass_Cl35
    non_halogen_mass = mass_molecular_ion - mass_of_two_Cl
    
    print("Step 1: Determine the mass of the hydrocarbon part of the molecule.")
    print(f"The molecular ion peak (M+) is at m/z = {mass_molecular_ion}.")
    print(f"This peak corresponds to the molecule with the lightest isotopes (two ³⁵Cl atoms).")
    print(f"Mass of two ³⁵Cl atoms = 2 * {mass_Cl35} = {mass_of_two_Cl} u.")
    print(f"Mass of the rest of the molecule = M+ - (mass of two ³⁵Cl atoms)")
    print(f"                             = {mass_molecular_ion} - {mass_of_two_Cl} = {non_halogen_mass} u.")
    
    # Propose the hydrocarbon formula
    num_C = 15
    num_H = 12
    calculated_hydrocarbon_mass = num_C * mass_C + num_H * mass_H
    
    print(f"\nA plausible formula for a mass of {non_halogen_mass} is C{num_C}H{num_H}, as {num_C}*{mass_C} + {num_H}*{mass_H} = {calculated_hydrocarbon_mass}.")
    print("This gives a total molecular formula of C₁₅H₁₂Cl₂.\n")

    # Step 2: Analyze the fragmentation pattern
    print("Step 2: Analyze the fragmentation pattern to confirm the structure.")
    # Loss of one chlorine atom
    calculated_fragment1 = mass_molecular_ion - mass_Cl35
    print(f"A major fragment is seen at m/z = {mass_M_minus_Cl}. This corresponds to the loss of one ³⁵Cl atom from the molecular ion:")
    print(f"[M - Cl]⁺ = {mass_molecular_ion} - {mass_Cl35} = {calculated_fragment1}.")
    
    # Formation of the base peak
    calculated_fragment2 = mass_M_minus_Cl - 2 * mass_H
    print(f"\nThe base peak (most intense) is at m/z = {mass_base_peak}. This is formed by the loss of a dihydrogen molecule (H₂) from the [M - Cl]⁺ fragment:")
    print(f"[M - Cl - H₂]⁺ = {mass_M_minus_Cl} - 2*{mass_H} = {calculated_fragment2}.")
    print("The high intensity of this peak suggests a very stable, rearranged cationic structure.\n")
    
    # Final Conclusion
    final_compound_name = "1,1-bis(4-chlorophenyl)prop-1-ene"
    print("Based on the molecular formula C₁₅H₁₂Cl₂ and the specific fragmentation pattern, the compound is identified.")
    
    print(f"\nFinal Answer: The IUPAC name of the compound is {final_compound_name}.")

solve_mass_spectrometry()