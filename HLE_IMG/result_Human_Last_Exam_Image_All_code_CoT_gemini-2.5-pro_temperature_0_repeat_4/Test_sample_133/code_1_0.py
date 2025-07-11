import re

def solve_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the compound.
    """
    # Data from the table in the image
    spectrum_data = {
        262: 29.0, 264: 12.9, 266: 3.22,
        225: 100.0, 227: 66.7,
        192: 24.5, 190: 48.4
    }
    
    # --- Step 1: Determine the Molecular Weight ---
    print("Step 1: Analyzing the Molecular Ion Peak (M+)")
    # The highest m/z cluster is at 262, 264, 266.
    molecular_ion_mass = 262
    print(f"The peak at the highest mass-to-charge ratio (m/z) is at {molecular_ion_mass}.")
    print(f"This suggests the molecular weight of the compound (containing the most common isotopes) is {molecular_ion_mass} u.\n")

    # --- Step 2: Analyze the Isotope Pattern ---
    print("Step 2: Analyzing the Isotope Pattern")
    m_peak_intensity = spectrum_data[262]
    m_plus_2_intensity = spectrum_data[264]
    m_plus_4_intensity = spectrum_data[266]
    
    # Normalize intensities relative to the M+ peak at 262
    relative_m_plus_2 = (m_plus_2_intensity / m_peak_intensity) * 100
    relative_m_plus_4 = (m_plus_4_intensity / m_peak_intensity) * 100
    
    # Theoretical ratio for 2 Chlorine atoms (35Cl:37Cl is approx 3:1) is M:M+2:M+4 ≈ 9:6:1 or 100:65:10
    print(f"The molecular ion cluster shows peaks at m/z {molecular_ion_mass}, {molecular_ion_mass+2}, and {molecular_ion_mass+4}.")
    print(f"Their intensity ratio is approximately 100 : {relative_m_plus_2:.0f} : {relative_m_plus_4:.0f}.")
    print("This pattern is characteristic of a compound containing two chlorine atoms.\n")
    num_cl = 2
    mass_cl = 35 # Using the mass of the most common isotope, 35Cl

    # --- Step 3: Deduce the Molecular Formula ---
    print("Step 3: Deducing the Molecular Formula")
    mass_of_cl_atoms = num_cl * mass_cl
    mass_of_rest = molecular_ion_mass - mass_of_cl_atoms
    print(f"The mass of the two chlorine atoms is 2 * {mass_cl} = {mass_of_cl_atoms} u.")
    print(f"The mass of the rest of the molecule (C_x H_y) is {molecular_ion_mass} - {mass_of_cl_atoms} = {mass_of_rest} u.")
    
    # Determine the formula for the C_x H_y part
    num_c = 15
    num_h = 12
    mass_c = 12
    mass_h = 1
    calculated_mass = num_c * mass_c + num_h * mass_h
    print(f"A plausible formula for a mass of {mass_of_rest} is C{num_c}H{num_h} (15 * 12 + 12 * 1 = {calculated_mass}).")
    print(f"Therefore, the proposed molecular formula is C{num_c}H{num_h}Cl{num_cl}.\n")

    # --- Step 4: Analyze the Fragmentation Pattern ---
    print("Step 4: Analyzing the Fragmentation Pattern")
    base_peak_mass = 225
    fragment_227_mass = 227
    fragment_192_mass = 192
    fragment_190_mass = 190

    print(f"The base peak (most intense fragment) is at m/z = {base_peak_mass}.")
    
    # Loss of Cl
    loss_cl = molecular_ion_mass - fragment_227_mass
    print(f"A major fragment at m/z {fragment_227_mass} corresponds to the loss of a fragment of mass {loss_cl} ({molecular_ion_mass} - {fragment_227_mass}). This is a chlorine atom (Cl).")
    print(f"Equation: [C{num_c}H{num_h}Cl{num_cl}]⁺˙ -> [C{num_c}H{num_h}Cl]⁺ + Cl˙")
    
    # Loss of Cl and H2
    loss_cl_h2 = molecular_ion_mass - base_peak_mass
    print(f"The base peak at m/z {base_peak_mass} corresponds to the loss of a fragment of mass {loss_cl_h2} ({molecular_ion_mass} - {base_peak_mass}). This can be explained by the loss of a chlorine atom ({mass_cl}) followed by the loss of a hydrogen molecule (H₂, mass 2).")
    print(f"Equation: [C{num_c}H{num_h}Cl]⁺ -> [C{num_c}H{num_h-2}Cl]⁺ + H₂")

    # Loss of 2Cl
    loss_2cl = molecular_ion_mass - fragment_192_mass
    print(f"The fragment at m/z {fragment_192_mass} corresponds to the loss of {loss_2cl} u, which is two chlorine atoms (2 * {mass_cl}).")
    print(f"Equation: [C{num_c}H{num_h}Cl{num_cl}]⁺˙ -> [C{num_c}H{num_h}]⁺˙ + 2Cl˙")

    # Loss of 2Cl and H2
    loss_2cl_h2 = molecular_ion_mass - fragment_190_mass
    print(f"The fragment at m/z {fragment_190_mass} corresponds to the loss of {loss_2cl_h2} u, which is two chlorine atoms and a hydrogen molecule.")
    print(f"Equation: [C{num_c}H{num_h}]⁺˙ -> [C{num_c}H{num_h-2}]⁺˙ + H₂\n")

    # --- Step 5: Identify the Compound ---
    print("Step 5: Identifying the Compound and Final Answer")
    print("The molecular formula C15H12Cl2 and the fragmentation pattern (loss of Cl, followed by H2) are consistent with the structure of 1,1-bis(4-chlorophenyl)prop-1-ene.")
    final_answer = "1,1'-prop-1-ene-1,1-diylbis(4-chlorobenzene)"
    print(f"\nThe compound is: {final_answer}")


solve_mass_spectrum()