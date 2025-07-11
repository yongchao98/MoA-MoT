def analyze_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the compound.
    """
    # Atomic masses (integer mass for simplicity)
    mass_C = 12
    mass_H = 1
    mass_Cl = 35

    # --- Step 1: Analyze Molecular Ion ---
    molecular_ion_mass = 260
    num_Cl = 3
    print(f"Step 1: Analysis of the Molecular Ion")
    print(f"The molecular ion peak (M+) for the all-³⁵Cl isotopologue is at m/z = {molecular_ion_mass}.")
    print(f"The isotopic cluster pattern suggests the presence of {num_Cl} chlorine atoms.\n")

    # --- Step 2: Determine Molecular Formula ---
    mass_of_3_Cl = num_Cl * mass_Cl
    mass_of_remainder = molecular_ion_mass - mass_of_3_Cl
    
    # Determine formula for the remainder (mass 155)
    num_C_remainder = 12
    num_H_remainder = 11
    calculated_mass_remainder = num_C_remainder * mass_C + num_H_remainder * mass_H
    
    print(f"Step 2: Molecular Formula Determination")
    print(f"Mass of {num_Cl} chlorine atoms = {num_Cl} * {mass_Cl} = {mass_of_3_Cl} u")
    print(f"Mass of the non-halogen part of the molecule = {molecular_ion_mass} - {mass_of_3_Cl} = {mass_of_remainder} u")
    print(f"A plausible formula for a mass of {mass_of_remainder} u is C{num_C_remainder}H{num_H_remainder} (mass = {calculated_mass_remainder} u).")
    print(f"Therefore, the molecular formula is C{num_C_remainder}H{num_H_remainder}Cl{num_Cl}.\n")

    # --- Step 3: Analyze Fragmentation ---
    base_peak_mass = 225
    fragment_143_mass = 143
    
    loss_for_base_peak = molecular_ion_mass - base_peak_mass
    mass_CCl3_radical = 1 * mass_C + 3 * mass_Cl
    loss_for_fragment_143 = molecular_ion_mass - fragment_143_mass

    print(f"Step 3: Fragmentation Analysis")
    print(f"The base peak is at m/z = {base_peak_mass}.")
    print(f"This corresponds to a loss of {molecular_ion_mass} - {base_peak_mass} = {loss_for_base_peak} u from the molecular ion.")
    print(f"This loss of {loss_for_base_peak} u is a chlorine atom, forming the [M-Cl]⁺ fragment.\n")
    
    print(f"A significant peak at m/z = {fragment_143_mass} is observed.")
    print(f"This corresponds to a loss of {molecular_ion_mass} - {fragment_143_mass} = {loss_for_fragment_143} u from the molecular ion.")
    print(f"This loss matches the mass of a trichloromethyl radical (•CCl₃), which is {mass_CCl3_radical} u.")
    print(f"This strongly suggests the presence of a -CCl₃ group in the molecule.\n")

    # --- Step 4: Propose Structure and Name ---
    # The structure is C11H11-CCl3. A plausible, though not common, structure is
    # 6,6,6-trichloro-1-phenyl-1,3-hexadiene.
    final_name = "6,6,6-trichloro-1-phenyl-1,3-hexadiene"
    print(f"Step 4: Structure and IUPAC Name")
    print(f"The evidence points to a structure of the type C₁₁H₁₁-CCl₃.")
    print(f"A plausible structure that fits the molecular formula C₁₂H₁₁Cl₃ is {final_name}.")

analyze_mass_spectrum()