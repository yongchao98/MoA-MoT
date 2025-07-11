def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized glycan ion.
    
    The process involves:
    1. Calculating the mass of the native A2G2S2 glycan.
    2. Adjusting the mass for the amidation of two sialic acids.
    3. Adjusting the mass for the permethylation of all available sites.
    4. Adding the mass of a sodium ion to get the final m/z for the [M+Na]+ ion.
    """
    
    # Monoisotopic atomic masses (Da)
    atomic_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Step 1: Calculate the mass of the native glycan (A2G2S2)
    # Composition: 2 GlcNAc, 3 Man, 2 Gal, 2 Neu5Ac
    # We calculate the mass from the sum of its monosaccharide residues plus one water molecule for the reducing end.
    
    # Mass of monosaccharide residues (the part within the polymer chain)
    # GlcNAc residue (C8H13NO5)
    mass_glcnac_res = 8 * atomic_mass['C'] + 13 * atomic_mass['H'] + 1 * atomic_mass['N'] + 5 * atomic_mass['O']
    # Hexose residue (Mannose/Galactose, C6H10O5)
    mass_hex_res = 6 * atomic_mass['C'] + 10 * atomic_mass['H'] + 5 * atomic_mass['O']
    # Sialic acid (Neu5Ac) residue (C11H17NO8)
    mass_neu5ac_res = 11 * atomic_mass['C'] + 17 * atomic_mass['H'] + 1 * atomic_mass['N'] + 8 * atomic_mass['O']
    # Mass of water (H2O) for the reducing end
    mass_h2o = 2 * atomic_mass['H'] + 1 * atomic_mass['O']

    # Total mass of the native glycan
    # It has 2 GlcNAc, 5 Hexose (3 Man + 2 Gal), and 2 Neu5Ac
    mass_native_glycan = (2 * mass_glcnac_res) + (5 * mass_hex_res) + (2 * mass_neu5ac_res) + mass_h2o

    # Step 2: Calculate mass change from amidation
    # The reaction R-COOH + NH3 -> R-CONH2 + H2O means the glycan loses an OH group and gains an NH2 group.
    # Mass change = Mass(NH2) - Mass(OH) = (N + 2H) - (O + H) = N + H - O
    mass_change_per_amidation = atomic_mass['N'] + atomic_mass['H'] - atomic_mass['O']
    num_sialic_acids = 2
    total_mass_change_amidation = num_sialic_acids * mass_change_per_amidation
    
    mass_amidated_glycan = mass_native_glycan + total_mass_change_amidation

    # Step 3: Calculate mass change from permethylation
    # A methyl group (CH3) replaces a hydrogen (H) at each site.
    # Mass change per site = Mass(CH3) - Mass(H) = Mass(CH2) = C + 2H
    mass_per_methylation = atomic_mass['C'] + 2 * atomic_mass['H']
    
    # Count methylation sites on the amidated glycan:
    # - Core (Man3GlcNAc2): 12 -OH groups
    # - GlcNAc N-acetyls: 2 -NH groups
    # - Galactose (x2): 2 * 3 = 6 -OH groups
    # - Amidated Sialic Acid (x2): Each has 4 -OH and 1 -NH (from N-acetyl). The new amide -CONH2 is not methylated.
    #   So, 2 * (4 + 1) = 10 sites.
    # Total sites = 12 + 2 + 6 + 10 = 30 sites
    num_methylation_sites = 30
    total_mass_increase_methylation = num_methylation_sites * mass_per_methylation

    # Step 4: Calculate the mass of the final derivatized molecule (M)
    mass_final_molecule = mass_amidated_glycan + total_mass_increase_methylation

    # Step 5: Calculate the m/z of the singly-sodiated ion [M+Na]+
    mass_na = atomic_mass['Na']
    final_mz = mass_final_molecule + mass_na

    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("Since the chemical modifications are not linkage-specific, they will all have the same final mass.")
    print("\nCalculation Breakdown:")
    print(f"Mass of the final derivatized glycan (M): {mass_final_molecule:.4f} Da")
    print(f"Mass of the sodium adduct (Na+): {mass_na:.4f} Da")
    print("\nFinal Equation and Result:")
    print(f"{mass_final_molecule:.4f} (M) + {mass_na:.4f} (Na+) = {final_mz:.4f} m/z")
    
    print("\nTherefore, the mass observed for all three glycans should be:")
    print(f"{final_mz:.4f}")
    
    return final_mz

# Execute the calculation and store the final answer
final_answer = calculate_glycan_mass()
# The final answer is returned in the format requested by the user.
# The value is rounded to 4 decimal places as is common in high-resolution mass spectrometry.
print(f'<<<{final_answer:.4f}>>>')
