def calculate_glycan_mass():
    """
    Calculates the theoretical m/z of a permethylated, amidated, and sodiated A2G2S2 glycan.
    """
    # Monoisotopic atomic masses
    MASS_H = 1.007825
    MASS_C = 12.000000
    MASS_O = 15.994915
    MASS_N = 14.003074
    MASS_Na = 22.989770

    # Step 1: Calculate the initial mass of the A2G2S2 glycan
    # Formula for A2G2S2 (HexNAc(4)Hex(5)NeuAc(2)): C84 H138 N6 O62
    num_C = 84
    num_H = 138
    num_N = 6
    num_O = 62
    initial_glycan_mass = (num_C * MASS_C) + (num_H * MASS_H) + (num_N * MASS_N) + (num_O * MASS_O)

    # Step 2: Calculate mass change from amidation
    # Two sialic acids, so two reactions.
    # Reaction: -COOH -> -CONH2. Net change per reaction: -O + N + H
    num_amidations = 2
    amidation_mass_change = num_amidations * (-MASS_O + MASS_N + MASS_H)

    # Step 3: Calculate mass change from permethylation
    # All free -OH and -NH groups are methylated (-H is replaced by -CH3).
    # Net change per methylation: +CH2
    # Number of sites: 32 (-OH) + 6 (N-acetyl -NH) = 38 on the native glycan.
    # Amidation removes one -OH and adds one -NH2 (two reactive H's) per sialic acid.
    # Net site change = -1 + 2 = +1 site per amidation.
    # Total sites = 38 + (2 * 1) = 40 sites.
    num_methylations = 40
    permethylation_mass_change = num_methylations * (MASS_C + (2 * MASS_H))

    # Step 4: Account for sodiation for [M+Na]+ ion
    sodiation_mass_change = MASS_Na
    
    # Step 5: Calculate the final m/z value
    final_mz = initial_glycan_mass + amidation_mass_change + permethylation_mass_change + sodiation_mass_change

    print("The three specified glycans (A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2) are isomers.")
    print("Isomers have the same chemical formula and therefore the same mass. The chemical modifications")
    print("also result in the same final mass. They will be separated by chromatography but will appear")
    print("at the same m/z in the mass spectrometer.\n")
    print("The calculation for the expected m/z is as follows:\n")
    
    print(f"1. Mass of initial glycan (C{num_C}H{num_H}N{num_N}O{num_O}): {initial_glycan_mass:.4f} Da")
    print(f"2. Mass change from amidation of {num_amidations} sialic acids: {amidation_mass_change:.4f} Da")
    print(f"3. Mass change from {num_methylations} permethylations: +{permethylation_mass_change:.4f} Da")
    print(f"4. Mass addition for sodiated ion [M+Na]+: +{sodiation_mass_change:.4f} Da")
    print("-" * 50)
    print("Final m/z Calculation:")
    print(f"{initial_glycan_mass:.4f} + ({amidation_mass_change:.4f}) + {permethylation_mass_change:.4f} + {sodiation_mass_change:.4f} = {final_mz:.4f}")

    # Return final value in the required format
    return f"<<<{final_mz:.4f}>>>"

# Run the calculation and print the final answer.
final_answer = calculate_glycan_mass()
print(f"\nThe single mass you should observe for all three isomeric glycans is {final_answer.strip('<>')}.")
print(final_answer)
