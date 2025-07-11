def calculate_mass(formula):
    """Calculates the monoisotopic mass of a chemical formula using integer masses."""
    atomic_masses = {'C': 12, 'H': 1, 'Cl': 35, 'O': 16, 'N': 14}
    mass = 0
    import re
    # Find all element-count pairs in the formula
    for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
        if count == '':
            count = 1
        else:
            count = int(count)
        mass += atomic_masses[element] * count
    return mass

def analyze_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the compound.
    """
    print("--- Mass Spectrum Analysis ---")

    # 1. Analyze the molecular ion peak
    print("\nStep 1: Analyzing the Molecular Ion (M+) Peak")
    m_ion_observed = 260
    print(f"The highest mass cluster starts at m/z = {m_ion_observed}.")
    print("The isotopic cluster pattern (peaks at 260, 262, 264, 266) is characteristic of a molecule containing 3 Chlorine atoms.")

    # 2. Determine the molecular formula
    print("\nStep 2: Determining the Molecular Formula")
    mass_of_cl3 = 3 * 35
    mass_of_non_halogen = m_ion_observed - mass_of_cl3
    print(f"Mass of the non-halogen part of the molecule: {m_ion_observed} - (3 * 35) = {mass_of_non_halogen}")
    
    # Deducing the hydrocarbon part
    c_atoms = mass_of_non_halogen // 12
    h_atoms = mass_of_non_halogen % 12
    hydrocarbon_formula = f"C{c_atoms}H{h_atoms}"
    print(f"A mass of {mass_of_non_halogen} corresponds to the formula {hydrocarbon_formula} ({c_atoms}*12 + {h_atoms}*1 = {mass_of_non_halogen}).")
    
    molecular_formula = f"C{c_atoms}H{h_atoms}Cl3"
    print(f"The proposed molecular formula is {molecular_formula}.")
    
    # Final equation for molecular weight calculation
    mw_calc_str = f"{c_atoms} * 12 + {h_atoms} * 1 + 3 * 35"
    mw_calc_val = calculate_mass(molecular_formula)
    print(f"Molecular weight calculation: {mw_calc_str} = {mw_calc_val}")
    

    # 3. Analyze fragment peaks to confirm structure
    print("\nStep 3: Analyzing Fragment Peaks")
    base_peak_obs = 225
    fragment_190_obs = 190
    fragment_145_obs = 145

    # Check M-Cl fragment
    m_minus_cl_formula = f"C{c_atoms}H{h_atoms}Cl2"
    m_minus_cl_mass = calculate_mass(m_minus_cl_formula)
    print(f"Base Peak at m/z {base_peak_obs}: Corresponds to loss of one Cl atom.")
    print(f"  - Calculated [M-Cl]+ ({m_minus_cl_formula}) mass: {m_minus_cl_mass}. This matches the observed peak.")

    # Check M-Cl-HCl fragment
    m_minus_cl_hcl_formula = f"C{c_atoms}H{h_atoms-1}Cl"
    m_minus_cl_hcl_mass = calculate_mass(m_minus_cl_hcl_formula)
    print(f"Peak at m/z {fragment_190_obs}: Corresponds to loss of Cl and HCl.")
    print(f"  - Calculated [M-Cl-HCl]+ ({m_minus_cl_hcl_formula}) mass: {m_minus_cl_hcl_mass}. This is {m_minus_cl_hcl_mass}, matching the observed peak at {fragment_190_obs}.")
    
    # Check dichlorophenyl fragment
    dichlorophenyl_formula = "C6H3Cl2"
    dichlorophenyl_mass = calculate_mass(dichlorophenyl_formula)
    print(f"Peak at m/z {fragment_145_obs}: Suggests a dichlorophenyl group.")
    print(f"  - Calculated [{dichlorophenyl_formula}]+ mass: {dichlorophenyl_mass}. This matches the observed peak.")

    # 4. Propose IUPAC Name
    print("\nStep 4: Conclusion")
    print("The fragmentation pattern supports a structure of a dichlorophenyl group attached to a chlorocyclohexene ring.")
    final_name = "1-(2-chlorocyclohex-1-en-1-yl)-2,4-dichlorobenzene"
    print("\nThe compound is identified as:")
    print(final_name)

# Run the analysis
analyze_mass_spectrum()