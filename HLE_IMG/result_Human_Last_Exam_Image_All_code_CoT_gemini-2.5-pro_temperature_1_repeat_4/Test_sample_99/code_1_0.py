def calculate_molecular_formula():
    """
    Calculates the molecular formula of the final product based on the reaction scheme.
    The final product is N-benzyl-2-((thiazol-2-yl)amino)-3-oxobutanamide.
    """

    # Atom counts for each part of the molecule
    # 1. Thiazole ring part: C3H3NS (thiazol-2-yl group)
    thiazole = {'C': 3, 'H': 2, 'N': 2, 'S': 1, 'O': 0} # The two ring CHs. The Ns are in the ring.

    # 2. Amino linker: -NH-
    amino_linker = {'C': 0, 'H': 1, 'N': 0, 'S': 0, 'O': 0} # The H on the amino group linking thiazole and the chain.

    # 3. Main chain: -CH(NH-thiazole)-C(=O)-CH3 from the butanoate part
    # And the new amide carbonyl: -C(=O)-NH-
    # Let's combine: CH3-C(=O)-CH-C(=O)-NH-
    main_chain = {'C': 4, 'H': 5, 'N': 1, 'S': 0, 'O': 2} # H count: CH3 (3) + CH (1) + NH (1) = 5

    # 4. Benzyl group: -CH2-C6H5
    benzyl_group = {'C': 7, 'H': 7, 'N': 0, 'S': 0, 'O': 0} # H count: CH2 (2) + C6H5 (5) = 7

    # Total atom counts
    C = thiazole['C'] + amino_linker['C'] + main_chain['C'] + benzyl_group['C']
    H = thiazole['H'] + amino_linker['H'] + main_chain['H'] + benzyl_group['H']
    N = thiazole['N'] + amino_linker['N'] + main_chain['N'] + benzyl_group['N']
    O = thiazole['O'] + amino_linker['O'] + main_chain['O'] + benzyl_group['O']
    S = thiazole['S'] + amino_linker['S'] + main_chain['S'] + benzyl_group['S']
    
    # Let's do a more structured calculation based on the full structure:
    # CH3-C(=O)-CH(NH-C3H2NS)-C(=O)-NH-CH2-C6H5
    # C atoms: 1(CH3) + 1(CO) + 1(CH) + 3(thiazole) + 1(CO) + 1(CH2) + 6(phenyl) = 14
    C_calc = 1 + 1 + 1 + 3 + 1 + 1 + 6
    # H atoms: 3(CH3) + 1(CH) + 1(NH) + 2(thiazole CH) + 1(amide NH) + 2(CH2) + 5(phenyl) = 15
    H_calc = 3 + 1 + 1 + 2 + 1 + 2 + 5
    # N atoms: 2(thiazole) + 1(amino linker) + 1(amide) = 3 (Note: the amino linker N is now counted in the amide group, let's recount)
    # N atoms: 2(thiazole ring) + 1(amide N) = 3
    N_calc = 2 + 1
    # O atoms: 1(ketone CO) + 1(amide CO) = 2
    O_calc = 1 + 1
    # S atoms: 1(thiazole)
    S_calc = 1
    
    print("Calculating the molecular formula of the final product, N-benzyl-2-((thiazol-2-yl)amino)-3-oxobutanamide:")
    print(f"Number of Carbon (C) atoms = 1 (from CH3) + 1 (from ketone C=O) + 1 (from CH) + 3 (from thiazole ring) + 1 (from amide C=O) + 1 (from CH2) + 6 (from phenyl ring) = {C_calc}")
    print(f"Number of Hydrogen (H) atoms = 3 (from CH3) + 1 (from CH) + 1 (from amino NH) + 2 (from thiazole ring) + 1 (from amide NH) + 2 (from CH2) + 5 (from phenyl ring) = {H_calc}")
    print(f"Number of Nitrogen (N) atoms = 2 (from thiazole ring) + 1 (from amide) = {N_calc}")
    print(f"Number of Oxygen (O) atoms = 1 (from ketone C=O) + 1 (from amide C=O) = {O_calc}")
    print(f"Number of Sulfur (S) atoms = 1 (from thiazole ring) = {S_calc}")

    molecular_formula = f"C{C_calc}H{H_calc}N{N_calc}O{O_calc}S{S_calc}"
    print(f"\nThe final molecular formula is: {molecular_formula}")

calculate_molecular_formula()