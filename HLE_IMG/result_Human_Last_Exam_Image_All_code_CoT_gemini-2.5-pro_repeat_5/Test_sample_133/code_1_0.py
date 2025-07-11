import math

def analyze_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the unknown compound.
    This function will print the step-by-step reasoning.
    """
    
    # --- Step 1 & 2: Analyze Key Peaks and Isotope Patterns ---
    print("Step 1: Analyzing the Base Peak and Molecular Ion")
    
    # Data from the table for key peaks
    base_peak = {'mz': 225, 'intensity': 100}
    base_peak_m_plus_2 = {'mz': 227, 'intensity': 66.7}
    
    molecular_ion_cluster = {
        258: 17.4,
        260: 34.8,
        262: 29.0,
        264: 12.9,
        266: 3.22
    }
    
    # Theoretical isotope ratios for Chlorine
    # For 2 Cl atoms, M:(M+2) ratio is approx. 100:64
    # For 3 Cl atoms, M:(M+2):(M+4) ratio is approx. 100:96:31
    
    # Analysis of the base peak at m/z 225
    ratio_base_peak = base_peak_m_plus_2['intensity'] / base_peak['intensity']
    print(f"The base peak is at m/z = {base_peak['mz']}.")
    print(f"The peak at m/z = {base_peak_m_plus_2['mz']} has an intensity of {base_peak_m_plus_2['intensity']}%.")
    print(f"The intensity ratio I(227)/I(225) is {ratio_base_peak:.2f}.")
    print("This is very close to the theoretical ratio of ~0.64 for a fragment containing two Chlorine atoms.")
    print("Conclusion: The fragment at m/z 225 contains 2 Chlorine atoms.\n")
    
    # Analysis of the molecular ion cluster
    # Let's assume the molecular ion M+ is the first peak of the highest-mass cluster, m/z 260.
    m_ion_mz = 260
    m_plus_2_mz = 262
    m_plus_4_mz = 264
    
    m_ion_intensity = molecular_ion_cluster[m_ion_mz]
    m_plus_2_intensity = molecular_ion_cluster[m_plus_2_mz]
    m_plus_4_intensity = molecular_ion_cluster[m_plus_4_mz]

    # Normalize to the M+ peak at 260
    normalized_m_plus_2 = (m_plus_2_intensity / m_ion_intensity) * 100
    normalized_m_plus_4 = (m_plus_4_intensity / m_ion_intensity) * 100

    print(f"The highest significant mass cluster starts at m/z = {m_ion_mz}.")
    print(f"The observed isotope pattern for this cluster (normalized to {m_ion_mz}=100) is roughly:")
    print(f"I({m_ion_mz}):100, I({m_plus_2_mz}):{normalized_m_plus_2:.0f}, I({m_plus_4_mz}):{normalized_m_plus_4:.0f}")
    print("This pattern (100:83:37) is a reasonable match for a molecule containing three Chlorine atoms (theoretical ~100:96:31).")
    print("Conclusion: The molecular ion is at m/z 260 and the molecule contains 3 Chlorine atoms.\n")
    
    # --- Step 3: Determine Fragmentation ---
    print("Step 3: Determining the Fragmentation Pattern")
    mass_loss = m_ion_mz - base_peak['mz']
    print(f"The mass difference between the molecular ion ({m_ion_mz}) and the base peak ({base_peak['mz']}) is {mass_loss} Da.")
    print("This mass loss corresponds to a single Chlorine atom (mass ~35 Da).")
    print("This confirms the fragmentation pathway: M⁺˙ -> [M-Cl]⁺.")
    print("This is consistent with a molecule containing 3 Cl atoms fragmenting to a base peak containing 2 Cl atoms.\n")
    
    # --- Step 4: Deduce the Molecular Formula ---
    print("Step 4: Deducing the Molecular Formula")
    mass_of_3_cl = 3 * 35
    mass_of_backbone = m_ion_mz - mass_of_3_cl
    print(f"The molecular weight of the all-³⁵Cl isotopologue is {m_ion_mz}.")
    print(f"Subtracting the mass of 3 Chlorine atoms (3 * 35 = {mass_of_3_cl}), the mass of the non-halogen part is {mass_of_backbone} Da.")
    print(f"A hydrocarbon formula for a mass of {mass_of_backbone} is C₁₂H₁₁ (12*12 + 11*1 = 144 + 11 = {mass_of_backbone}).")
    print("Therefore, the molecular formula of the compound is C₁₂H₁₁Cl₃.\n")
    
    # --- Step 5: Identify the Compound ---
    print("Step 5: Identifying the Compound and its IUPAC Name")
    print("The formula C₁₂H₁₁Cl₃ has a degree of unsaturation of 6.")
    print("This suggests a structure containing one aromatic ring (DBE=4) and two additional double bonds or a second ring (DBE=2).")
    print("A plausible structure is a phenyl-substituted trichlorocyclohexadiene.")
    print("Based on the fragmentation pattern, a likely candidate is an isomer of phenyltrichlorocyclohexadiene.")
    
    final_answer = "4,5,6-trichloro-1-phenylcyclohexa-1,3-diene"
    print(f"\nThe compound is an isomer of phenyltrichlorocyclohexadiene. A specific IUPAC name for such a structure is:")
    print(final_answer)

analyze_mass_spectrum()