def solve_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the compound.
    """
    # Key data points from the provided spectrum table
    peaks = {
        266: 3.22, 264: 12.9, 262: 29.0, 260: 34.8,
        229: 22.2, 227: 66.7, 225: 100.0,
        192: 24.5, 190: 48.4,
        155: 16.1,
        147: 1.02, 143: 27.8
    }

    # --- Step 1 & 2: Analyze Molecular Ion (M+) and determine Cl count ---
    print("Step 1: Analyzing the Molecular Ion (M+) Region")
    m_ion_cluster_obs = {
        '260': peaks[260], '262': peaks[262], '264': peaks[264], '266': peaks[266]
    }
    # Normalize to the most intense peak in the cluster (m/z 260)
    norm_factor = m_ion_cluster_obs['260']
    m_ion_ratio_obs = {k: round(100 * v / norm_factor, 1) for k, v in m_ion_cluster_obs.items()}
    
    # Theoretical ratios for 3 Chlorine atoms (M, M+2, M+4, M+6)
    m_ion_ratio_theo_3cl = {'M': 100, 'M+2': 97.5, 'M+4': 32.5, 'M+6': 3.7}
    
    print(f"The highest mass cluster is observed at m/z = {list(m_ion_cluster_obs.keys())}.")
    print(f"Observed intensity ratio: M(260):M+2(262):M+4(264):M+6(266) ≈ {m_ion_ratio_obs['260']}:{m_ion_ratio_obs['262']}:{m_ion_ratio_obs['264']}:{m_ion_ratio_obs['266']}")
    print(f"Theoretical ratio for 3 Cl atoms: {m_ion_ratio_theo_3cl['M']}:{m_ion_ratio_theo_3cl['M+2']}:{m_ion_ratio_theo_3cl['M+4']}:{m_ion_ratio_theo_3cl['M+6']}")
    print("Conclusion: The pattern strongly suggests the presence of 3 Chlorine atoms and a molecular weight of 260 amu (for the all-³⁵Cl isotopomer).\n")

    # --- Step 3: Identify Base Peak and Fragmentation Cascade ---
    print("Step 2: Analyzing Fragmentation")
    molecular_ion_mass = 260
    base_peak_mass = 225
    mass_loss_cl = 35

    print(f"The base peak (highest intensity) is at m/z = {base_peak_mass}.")
    print(f"Mass difference from M+ ({molecular_ion_mass}) = {molecular_ion_mass} - {base_peak_mass} = {molecular_ion_mass - base_peak_mass} amu.")
    print(f"This corresponds to the loss of one Chlorine atom (mass ≈ {mass_loss_cl}). The base peak is [M-Cl]⁺.")
    
    # Check isotope pattern of base peak
    base_peak_m_plus_2 = 227
    ratio_base_peak = peaks[base_peak_m_plus_2] / peaks[base_peak_mass]
    print(f"The base peak's own isotope pattern, I(227)/I(225) = {peaks[base_peak_m_plus_2]}/{peaks[base_peak_mass]} = {ratio_base_peak:.2f}, matches the theoretical ratio for 2 Cl atoms (~0.65).")
    
    fragment_m_minus_2cl = 190
    fragment_m_minus_3cl = 155
    print("A clear fragmentation cascade of sequential chlorine loss is observed:")
    print(f"  M⁺ (m/z {molecular_ion_mass}) -> [M-Cl]⁺ (m/z {base_peak_mass}) -> [M-2Cl]⁺ (m/z {fragment_m_minus_2cl}) -> [M-3Cl]⁺ (m/z {fragment_m_minus_3cl})\n")

    # --- Step 4 & 5: Determine Formula and Structure ---
    print("Step 3: Determining Formula and Structure")
    core_fragment_mass = 155
    print(f"The chlorine-free core fragment [M-3Cl]⁺ has m/z = {core_fragment_mass}.")
    print("Based on fragmentation patterns (aliphatic-like loss of Cl, lack of aromatic fragments), the most likely formula for a mass of 260 with 3 Cl is C₁₀H₁₉OCl₃.")
    
    print("\nProposing a structure based on further fragmentation:")
    print("The structure 1,1,1-Trichlorodecan-2-ol is proposed.")
    alpha_cleavage_frag_1 = 143 # Loss of CCl3 radical (117) -> 260 - 117 = 143
    alpha_cleavage_frag_2 = 147 # Loss of C8H17 radical (113) -> 260 - 113 = 147
    print(f"This structure would undergo alpha-cleavage characteristic of a secondary alcohol:")
    print(f"  - Loss of ·CCl₃ would produce a fragment at m/z = {molecular_ion_mass} - 117 = {alpha_cleavage_frag_1}. A peak is observed at m/z {alpha_cleavage_frag_1}.")
    print(f"  - Loss of ·C₈H₁₇ would produce a fragment at m/z = {molecular_ion_mass} - 113 = {alpha_cleavage_frag_2}. A peak is observed at m/z {alpha_cleavage_frag_2}.")
    print("The presence of these fragments strongly supports this structural assignment.\n")
    
    # --- Final Answer ---
    final_answer = "1,1,1-Trichlorodecan-2-ol"
    print("Final Answer: The IUPAC name of the compound is:")
    print(final_answer)

solve_mass_spectrum()