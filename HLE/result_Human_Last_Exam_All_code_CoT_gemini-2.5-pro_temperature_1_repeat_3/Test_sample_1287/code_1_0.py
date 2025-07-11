def calculate_glycan_mass():
    """
    Calculates the m/z of a derivatized sialylated N-glycan.
    
    The steps are:
    1. Start with the mass of the native glycan (A2G2S2).
    2. Apply mass change from amidation of two sialic acids.
    3. Apply mass change from permethylation.
    4. Add the mass of a sodium ion for the final m/z.
    """

    # Monoisotopic masses of atoms and groups
    # Using high precision values for accuracy
    H = 1.007825
    C = 12.000000
    N = 14.003074
    O = 15.994915
    Na = 22.989770

    # Mass of groups involved in the reactions
    mass_OH = O + H
    mass_NH2 = N + 2 * H
    mass_CH2 = C + 2 * H # Net addition during methylation (CH3 replaces H)

    # 1. Starting mass of the native glycan (A2G2S2: Hex(5)HexNAc(4)NeuAc(2))
    # This mass is the same for all three of your isomeric glycans.
    mass_native_glycan = 2694.9458

    # 2. Amidation of two sialic acids (-COOH becomes -CONH2)
    # The mass change is the mass of NH2 minus the mass of OH.
    num_sialic_acids = 2
    mass_after_amidation = mass_native_glycan - num_sialic_acids * mass_OH + num_sialic_acids * mass_NH2

    # 3. Permethylation
    # We need to count the number of methylation sites (-OH and -NH groups).
    # The number of sites is constant regardless of the sialic acid linkage (2,3 vs 2,6).
    # Breakdown of sites on the amidated glycan:
    # - Reducing GlcNAc: 3 OH + 1 NH = 4
    # - Core GlcNAc: 2 OH + 1 NH = 3
    # - Branching Man: 2 OH = 2
    # - Arm Man (x2): 2 * 3 OH = 6
    # - Arm GlcNAc (x2): 2 * (2 OH + 1 NH) = 6
    # - Arm Gal (x2): 2 * 3 OH = 6
    # - Terminal Sia-Amide (x2): 2 * (4 OH + 1 NH(acetyl) + 1 CONH2(2 sites)) = 14
    # Total sites = 4 + 3 + 2 + 6 + 6 + 6 + 14 = 41
    num_methyl_sites = 41
    mass_increase_from_methylation = num_methyl_sites * mass_CH2
    mass_after_permethylation = mass_after_amidation + mass_increase_from_methylation

    # 4. Formation of singly-sodiated ion ([M+Na]+)
    final_mass_ion = mass_after_permethylation + Na
    
    # 5. Final m/z for a +1 ion
    final_mz = final_mass_ion / 1

    # Print the detailed calculation
    print("Since the chemical modifications are not specific to the sialic acid linkage, all three of your glycans are expected to have the same mass-to-charge ratio.")
    print("\nThe calculation is as follows:\n")
    print("Equation: Mass(Native) - 2*Mass(OH) + 2*Mass(NH2) + 41*Mass(CH2) + Mass(Na)\n")
    print("Calculation:")
    print(f"{mass_native_glycan:.4f} - 2*{mass_OH:.4f} + 2*{mass_NH2:.4f} + {num_methyl_sites}*{mass_CH2:.4f} + {Na:.4f} = {final_mz:.4f}")
    
    print("\nThe expected m/z for A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 is:")
    print(f"{final_mz:.4f}")
    
    return final_mz

# Execute the calculation
final_mz_value = calculate_glycan_mass()
# The final answer format is requested at the very end
# <<<3290.6092>>>