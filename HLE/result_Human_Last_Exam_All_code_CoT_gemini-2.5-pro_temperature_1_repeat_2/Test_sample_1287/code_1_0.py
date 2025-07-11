def calculate_glycan_masses():
    """
    Calculates the expected m/z for three differentially sialylated and derivatized glycans.
    """
    # Step 1: Define monoisotopic masses of atoms and molecules
    mass_C = 12.000000
    mass_H = 1.007825
    mass_O = 15.994915
    mass_N = 14.003074
    mass_Na = 22.989770

    # Step 2: Calculate the mass of the underivatized base glycan: A2G2S2
    # Composition: 4x GlcNAc (C8H15NO6), 5x Hexose (C6H12O6), 2x Neu5Ac (C11H19NO9)
    # To get the polymer mass, we sum the monomers and subtract water (H2O) for each of the 10 glycosidic bonds.
    mass_glcnac = 8*mass_C + 15*mass_H + 1*mass_N + 6*mass_O
    mass_hex = 6*mass_C + 12*mass_H + 6*mass_O
    mass_neuac = 11*mass_C + 19*mass_H + 1*mass_N + 9*mass_O
    mass_h2o = 2*mass_H + 1*mass_O
    
    # Total 11 residues means 10 bonds for a biantennary glycan.
    mass_underivatized_glycan = (4 * mass_glcnac + 5 * mass_hex + 2 * mass_neuac) - (10 * mass_h2o)

    # Step 3: Calculate the mass for the reference glycan, A2G(4)2S(3)2
    # Both sialic acids are a-2,3, so they are not amidated but are methyl-esterified.
    # The entire glycan is permethylated.
    # A standard A2G2S2 glycan has 32 sites for methylation (OH, NH, and the two COOH groups).
    # Each methylation replaces H with CH3, so the mass increase is CH2.
    num_methylations = 32
    mass_ch2 = mass_C + 2 * mass_H
    mass_increase_from_methylation = num_methylations * mass_ch2
    
    mass_glycan1_neutral = mass_underivatized_glycan + mass_increase_from_methylation
    mass_glycan1_sodiated = mass_glycan1_neutral + mass_Na
    
    print("Calculation for A2G(4)2S(3)2 (both sialic acids are methyl-esterified):")
    print(f"Base Glycan Mass + Methylation Mass + Sodium Mass = m/z")
    print(f"{mass_underivatized_glycan:.3f} + {mass_increase_from_methylation:.3f} + {mass_Na:.3f} = {mass_glycan1_sodiated:.3f}")
    print("-" * 20)

    # Step 4: Calculate the mass difference (delta) for the other two glycans.
    # The difference is replacing a methyl esterification with an amidation.
    # Methyl esterification of -COOH adds a CH2 group (net change from -COOH to -COOCH3 replaces H with CH3).
    # Amidation of -COOH to -CONH2 replaces an O with an NH.
    mass_change_esterification = mass_ch2
    mass_change_amidation = mass_N + mass_H - mass_O # Replaces OH with NH2
    
    delta_amide_vs_ester = mass_change_amidation - mass_change_esterification
    
    # Step 5: Calculate the mass for A2G(4)S(3)S(6)
    # One sialic acid is amidated, so we apply the delta once.
    mass_glycan2_sodiated = mass_glycan1_sodiated + delta_amide_vs_ester
    
    print("Calculation for A2G(4)S(3)S(6) (one ester, one amide):")
    print(f"Reference m/z + Mass Difference = m/z")
    print(f"{mass_glycan1_sodiated:.3f} + ({delta_amide_vs_ester:.3f}) = {mass_glycan2_sodiated:.3f}")
    print("-" * 20)

    # Step 6: Calculate the mass for A2G(4)2S(6)2
    # Both sialic acids are amidated, so we apply the delta twice.
    mass_glycan3_sodiated = mass_glycan1_sodiated + 2 * delta_amide_vs_ester
    
    print("Calculation for A2G(4)2S(6)2 (both sialic acids are amidated):")
    print(f"Reference m/z + 2 * (Mass Difference) = m/z")
    print(f"{mass_glycan1_sodiated:.3f} + 2 * ({delta_amide_vs_ester:.3f}) = {mass_glycan3_sodiated:.3f}")

calculate_glycan_masses()