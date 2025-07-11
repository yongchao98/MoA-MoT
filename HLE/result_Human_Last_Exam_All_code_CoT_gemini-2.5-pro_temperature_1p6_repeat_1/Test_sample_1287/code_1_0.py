def calculate_glycan_mass():
    """
    Calculates the expected m/z for specific sialylated glycans after
    amidation, permethylation, and sodiation.
    """
    # Monoisotopic atomic masses
    C = 12.000000
    H = 1.007825
    O = 15.994915
    N = 14.003074
    Na = 22.989770

    # Step 1: Define the glycan composition and calculate its initial mass.
    # The structure is Neu5Ac(2)Gal(2)Man(3)GlcNAc(4).
    # Molecular formulas for the full monosaccharides:
    # Neu5Ac: C11 H19 N O9
    # Gal: C6 H12 O6
    # Man: C6 H12 O6
    # GlcNAc: C8 H15 N O6

    num_neu5ac = 2
    num_gal = 2
    num_man = 3
    num_glcnac = 4

    total_monosaccharides = num_neu5ac + num_gal + num_man + num_glcnac
    num_glycosidic_bonds = total_monosaccharides - 1

    # Calculate the formula of the intact glycan (sum of all atoms minus water for each bond)
    base_formula_C = (num_neu5ac * 11) + (num_gal * 6) + (num_man * 6) + (num_glcnac * 8)
    base_formula_H = (num_neu5ac * 19) + (num_gal * 12) + (num_man * 12) + (num_glcnac * 15) - (num_glycosidic_bonds * 2)
    base_formula_N = (num_neu5ac * 1) + (num_glcnac * 4)
    base_formula_O = (num_neu5ac * 9) + (num_gal * 6) + (num_man * 6) + (num_glcnac * 6) - num_glycosidic_bonds

    initial_mass = (base_formula_C * C) + (base_formula_H * H) + (base_formula_N * N) + (base_formula_O * O)

    print("--- Step 1: Initial Glycan Mass Calculation ---")
    print(f"The glycan composition is Neu5Ac({num_neu5ac})Gal({num_gal})Man({num_man})GlcNAc({num_glcnac}).")
    print(f"This corresponds to the molecular formula: C{base_formula_C}H{base_formula_H}N{base_formula_N}O{base_formula_O}")
    print(f"Initial Mass = {initial_mass:.6f} Da\n")

    # Step 2: Account for the amidation of the two sialic acids.
    # The reaction converts -COOH to -CONH2. Net change is -OH + NH2 = -O +NH.
    mass_change_per_amidation = (N + H) - O
    total_amidation_change = num_neu5ac * mass_change_per_amidation
    mass_after_amidation = initial_mass + total_amidation_change

    print("--- Step 2: Amidation Mass Calculation ---")
    print(f"Amidation converts 2 -COOH groups to -CONH2 groups.")
    print(f"Mass change per amidation (-O +NH) = {N:.6f} + {H:.6f} - {O:.6f} = {mass_change_per_amidation:.6f} Da")
    print(f"Mass after amidation = {initial_mass:.6f} + 2 * ({mass_change_per_amidation:.6f}) = {mass_after_amidation:.6f} Da\n")

    # Step 3: Account for permethylation.
    # The key point is that changing the sialic acid linkage from alpha-2,3 to alpha-2,6
    # just moves the attachment point on the galactose residue (from C3-OH to C6-OH).
    # It does not change the total number of free -OH groups on the galactose or the glycan as a whole.
    # All three isomers have the same number of methylation sites.
    # Based on the standard biantennary glycan structure, there are 36 sites for methylation:
    # 2xSia (5 sites each) + 2xGal (3 sites each) + 2xAntenna GlcNAc (3 sites each) +
    # Core (3+3+1 for Mans, 3+4 for GlcNAcs) = 10+6+6+7+7 = 36 sites.
    num_methylation_sites = 36
    # Permethylation replaces -H with -CH3. Net change is +CH2.
    mass_change_per_methylation = C + (2 * H)
    total_methylation_change = num_methylation_sites * mass_change_per_methylation
    mass_after_methylation = mass_after_amidation + total_methylation_change

    print("--- Step 3: Permethylation Mass Calculation ---")
    print("The sialic acid linkage (2,3 vs 2,6) does not change the total number of available -OH and -NH groups.")
    print(f"All three isomers have {num_methylation_sites} methylation sites.")
    print(f"Mass change per methylation (+CH2) = {C:.6f} + 2*{H:.6f} = {mass_change_per_methylation:.6f} Da")
    print(f"Mass after permethylation = {mass_after_amidation:.6f} + {num_methylation_sites} * {mass_change_per_methylation:.6f} = {mass_after_methylation:.6f} Da\n")

    # Step 4: Account for sodiation for MS detection ([M+Na]+).
    final_mass = mass_after_methylation + Na
    
    print("--- Step 4: Sodiated Ion Mass Calculation ---")
    print(f"The mass spectrometer detects a singly-sodiated ion ([M+Na]+).")
    print(f"Final Observed m/z = {mass_after_methylation:.6f} (M) + {Na:.6f} (Na) = {final_mass:.6f} Da\n")
    
    # Final answer
    print("--- Final Answer ---")
    print("Because all three glycans are isomers and the chemical modifications result in the same mass change for each, their expected masses are identical.")
    print(f"The expected mass for A2G(4)2S(3)2 is: {final_mass:.4f}")
    print(f"The expected mass for A2G(4)S(3)S(6) is: {final_mass:.4f}")
    print(f"The expected mass for A2G(4)2S(6)2 is: {final_mass:.4f}")

    # Hidden final answer tag
    # return f"<<<{final_mass:.4f}, {final_mass:.4f}, {final_mass:.4f}>>>"


if __name__ == '__main__':
    calculate_glycan_mass()