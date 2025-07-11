def calculate_glycan_masses():
    """
    Calculates the expected m/z for three derivatized and permethylated glycans.
    """
    # Monoisotopic atomic masses
    MASS = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770
    }

    # Base composition of A2G2S2 glycan: (Neu5Ac)2(Gal)2(Man)3(GlcNAc)4
    # Formula determined by summing monomers and subtracting 10*H2O for linkages.
    # C(2*11+5*6+4*8) H(2*19+5*12+4*15-10*2) N(2*1+4*1) O(2*9+5*6+4*6-10*1)
    base_formula = {'C': 84, 'H': 138, 'N': 6, 'O': 62}
    
    # Baseline number of methylation sites for a standard permethylated A2G2S2,
    # where carboxyls are converted to methyl esters. This is a known value.
    baseline_methylation_sites = 33

    print("Analysis of Derivatized Sialylated Glycans")
    print("-" * 50)
    print(f"Starting with base A2G2S2 glycan with formula: C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}\n")

    # --- Glycan 1: A2G(4)2S(3)2 ---
    print("1. Glycan A2G(4)2S(3)2 (two α-2,3 linkages)")
    print("   Reaction: Both sialic acids undergo lactonization.")
    
    # Formula change for lactonization: -H2O per reaction
    formula_g1 = base_formula.copy()
    formula_g1['H'] -= 2 * 2  # 2x H2O loss
    formula_g1['O'] -= 2 * 1  # 2x H2O loss

    # Methylation sites: Lactonization consumes one COOH and one OH site.
    # For two lactones, 4 baseline sites are removed.
    methylations_g1 = baseline_methylation_sites - 4
    
    # Add methyl groups (add CH3, remove H for each site)
    formula_g1['C'] += methylations_g1 * 1
    formula_g1['H'] += methylations_g1 * 2

    mass_g1 = (formula_g1['C'] * MASS['C'] + 
               formula_g1['H'] * MASS['H'] + 
               formula_g1['N'] * MASS['N'] + 
               formula_g1['O'] * MASS['O'])
    mz_g1 = mass_g1 + MASS['Na']

    print(f"   Final Formula: C{formula_g1['C']}H{formula_g1['H']}N{formula_g1['N']}O{formula_g1['O']}")
    print(f"   Mass Equation: ({formula_g1['C']} * {MASS['C']:.6f}) + ({formula_g1['H']} * {MASS['H']:.6f}) + ({formula_g1['N']} * {MASS['N']:.6f}) + ({formula_g1['O']} * {MASS['O']:.6f}) + {MASS['Na']:.6f}")
    print(f"   Expected m/z for [M+Na]+: {mz_g1:.4f}\n")

    # --- Glycan 2: A2G(4)S(3)S(6) ---
    print("2. Glycan A2G(4)S(3)S(6) (one α-2,3 and one α-2,6 linkage)")
    print("   Reaction: One lactonization and one amidation.")
    
    # Formula change: 1x lactonization (-H2O) and 1x amidation (-O +NH)
    formula_g2 = base_formula.copy()
    formula_g2['H'] -= 2  # from lactonization
    formula_g2['O'] -= 1  # from lactonization
    formula_g2['O'] -= 1  # from amidation
    formula_g2['N'] += 1  # from amidation
    formula_g2['H'] += 1  # from amidation

    # Methylation sites: Lactone removes 2 sites. Amidation removes 1 (COOH) but adds 2 (NH2), a net gain of 1.
    # Total change = -2 (lactone) + 1 (amide) = -1 site.
    methylations_g2 = baseline_methylation_sites - 1
    
    # Add methyl groups
    formula_g2['C'] += methylations_g2 * 1
    formula_g2['H'] += methylations_g2 * 2
    
    mass_g2 = (formula_g2['C'] * MASS['C'] + 
               formula_g2['H'] * MASS['H'] + 
               formula_g2['N'] * MASS['N'] + 
               formula_g2['O'] * MASS['O'])
    mz_g2 = mass_g2 + MASS['Na']

    print(f"   Final Formula: C{formula_g2['C']}H{formula_g2['H']}N{formula_g2['N']}O{formula_g2['O']}")
    print(f"   Mass Equation: ({formula_g2['C']} * {MASS['C']:.6f}) + ({formula_g2['H']} * {MASS['H']:.6f}) + ({formula_g2['N']} * {MASS['N']:.6f}) + ({formula_g2['O']} * {MASS['O']:.6f}) + {MASS['Na']:.6f}")
    print(f"   Expected m/z for [M+Na]+: {mz_g2:.4f}\n")

    # --- Glycan 3: A2G(4)2S(6)2 ---
    print("3. Glycan A2G(4)2S(6)2 (two α-2,6 linkages)")
    print("   Reaction: Both sialic acids undergo amidation.")
    
    # Formula change for amidation: -O +NH per reaction
    formula_g3 = base_formula.copy()
    formula_g3['O'] -= 2 * 1 # 2x amidation
    formula_g3['N'] += 2 * 1 # 2x amidation
    formula_g3['H'] += 2 * 1 # 2x amidation
    
    # Methylation sites: Each amidation removes 1 COOH site but adds 2 NH2 sites (net +1 per amide).
    # Total change = +2 sites
    methylations_g3 = baseline_methylation_sites + 2

    # Add methyl groups
    formula_g3['C'] += methylations_g3 * 1
    formula_g3['H'] += methylations_g3 * 2

    mass_g3 = (formula_g3['C'] * MASS['C'] + 
               formula_g3['H'] * MASS['H'] + 
               formula_g3['N'] * MASS['N'] + 
               formula_g3['O'] * MASS['O'])
    mz_g3 = mass_g3 + MASS['Na']

    print(f"   Final Formula: C{formula_g3['C']}H{formula_g3['H']}N{formula_g3['N']}O{formula_g3['O']}")
    print(f"   Mass Equation: ({formula_g3['C']} * {MASS['C']:.6f}) + ({formula_g3['H']} * {MASS['H']:.6f}) + ({formula_g3['N']} * {MASS['N']:.6f}) + ({formula_g3['O']} * {MASS['O']:.6f}) + {MASS['Na']:.6f}")
    print(f"   Expected m/z for [M+Na]+: {mz_g3:.4f}\n")
    
    return f"{mz_g1:.2f}, {mz_g2:.2f}, {mz_g3:.2f}"


# Execute the calculation and format the final answer.
final_masses = calculate_glycan_masses()
# print(f"\nFinal Answer: {final_masses}")