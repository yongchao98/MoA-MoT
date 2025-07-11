def solve_mass_spectrometry_puzzle():
    """
    This script formalizes the step-by-step analysis of the mass spectrum to deduce the compound's identity.
    """
    # 1. Key data points from the spectrum
    molecular_ion_M = 262  # Corresponds to C15H12(35Cl)2
    molecular_ion_M_plus_2 = 264
    molecular_ion_M_plus_4 = 266
    
    base_peak = 225
    fragment_M_minus_Cl = 227 # Corresponds to [C15H12(35Cl)]+
    fragment_M_minus_2Cl_H2 = 190 # Corresponds to [C15H10]+
    
    # 2. Molecular formula deduction
    mass_of_hydrocarbon_part = molecular_ion_M - 2 * 35
    molecular_formula = "C15H12Cl2"
    
    # 3. Fragmentation pathway
    pathway_step_1 = f"{molecular_ion_M} (M+) -> {fragment_M_minus_Cl} ([M-Cl]+)"
    pathway_step_2 = f"{fragment_M_minus_Cl} ([M-Cl]+) -> {base_peak} ([M-Cl-H2]+, Base Peak)"
    pathway_step_3 = f"{base_peak} ([M-Cl-H2]+) -> {fragment_M_minus_2Cl_H2} ([M-2Cl-H2]+)"
    
    # 4. Final identification
    iupac_name = "1,1-dichloro-2,3-diphenylcyclopropane"

    print("Step-by-step deduction from the mass spectrum:")
    print("1. Molecular ion cluster found at m/z 262, 264, 266, suggesting a C(x)H(y)Cl2 compound.")
    print(f"2. Mass of the non-halogen part (C(x)H(y)) is {molecular_ion_M} - 70 = {mass_of_hydrocarbon_part} u.")
    print(f"3. A plausible formula for this mass is C15H12. This gives a molecular formula of {molecular_formula}.")
    print("4. The fragmentation pattern confirms this hypothesis:")
    print(f"   - Loss of Cl: {pathway_step_1}")
    print(f"   - Loss of H2 to form the stable base peak: {pathway_step_2}")
    print(f"   - Loss of the second Cl: {pathway_step_3}")
    print("\nConclusion:")
    print(f"The compound corresponding to the spectrum is {iupac_name}.")

solve_mass_spectrometry_puzzle()