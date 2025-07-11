def calculate_glycan_mass():
    """
    Calculates the m/z of a permethylated, amidated, and sodiated A2G2S2 glycan.
    """
    # Monoisotopic atomic masses
    ATOMIC_MASS = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770
    }

    # Monosaccharide compositions (as if they were free molecules)
    MONOSACCHARIDE_FORMULA = {
        'GlcNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6},
        'Man':    {'C': 6, 'H': 12, 'N': 0, 'O': 6},
        'Gal':    {'C': 6, 'H': 12, 'N': 0, 'O': 6},
        'Neu5Ac': {'C': 11, 'H': 19, 'N': 1, 'O': 9} # N-acetylneuraminic acid
    }
    
    # Water formula to subtract for glycosidic bonds
    H2O = {'C': 0, 'H': 2, 'N': 0, 'O': 1}

    # Glycan composition
    glycan_comp = {'GlcNAc': 4, 'Man': 3, 'Gal': 2, 'Neu5Ac': 2}
    
    # 1. Start with the sum of all monosaccharides
    base_formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for saccharide, count in glycan_comp.items():
        for atom, num in MONOSACCHARIDE_FORMULA[saccharide].items():
            base_formula[atom] += num * count
    
    # Subtract water for glycosidic linkages
    num_residues = sum(glycan_comp.values())
    num_links = num_residues - 1
    
    final_formula = base_formula.copy()
    for atom, num in H2O.items():
        final_formula[atom] -= num * num_links

    print("Step 1: Base Glycan (A2G2S2) Calculation")
    print(f"   - Formula of naked glycan polymer: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")

    # 2. Amidation of 2 sialic acids (-COOH -> -CONH2)
    # Net change per reaction is -O, +N, +H
    amidation_change = {'C': 0, 'H': 1 * 2, 'N': 1 * 2, 'O': -1 * 2}
    for atom, change in amidation_change.items():
        final_formula[atom] += change
        
    print("\nStep 2: Amidation")
    print("   - Two sialic acids are amidated (-COOH -> -CONH2).")
    print(f"   - Formula after amidation: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")

    # 3. Permethylation (replace H with CH3 on all OH and NH groups)
    # Net change per methylation is +CH2
    # The number of methylation sites is determined from structural analysis.
    # For a released, amidated A2G2S2 glycan, there are 41 sites.
    # This includes OH groups, acetyl-NH groups, and the new amide-NH2 groups.
    num_methylations = 41
    methylation_change = {'C': 1 * num_methylations, 'H': 2 * num_methylations, 'N': 0, 'O': 0}
    for atom, change in methylation_change.items():
        final_formula[atom] += change
        
    print("\nStep 3: Permethylation")
    print(f"   - Assumed number of methylation sites: {num_methylations}")
    print(f"   - Formula after permethylation: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")

    # 4. Calculate the mass of the derivatized neutral molecule
    neutral_mass = 0
    mass_eq_parts = []
    for atom, count in final_formula.items():
        mass = count * ATOMIC_MASS[atom]
        neutral_mass += mass
        mass_eq_parts.append(f"{count}*{ATOMIC_MASS[atom]}[{atom}]")

    print("\nStep 4: Calculate Mass of Neutral Derivatized Glycan [M]")
    print("   - Mass [M] = " + " + ".join(mass_eq_parts))
    print(f"   - Mass [M] = {final_formula['C']*ATOMIC_MASS['C']:.5f} + {final_formula['H']*ATOMIC_MASS['H']:.5f} + {final_formula['N']*ATOMIC_MASS['N']:.5f} + {final_formula['O']*ATOMIC_MASS['O']:.5f}")
    print(f"   - Calculated neutral mass [M]: {neutral_mass:.5f} Da")

    # 5. Add mass of sodium
    sodiated_mass = neutral_mass + ATOMIC_MASS['Na']

    print("\nStep 5: Add Sodium Ion for [M+Na]+")
    print(f"   - Final m/z [M+Na]+ = {neutral_mass:.5f} [M] + {ATOMIC_MASS['Na']} [Na]")
    print(f"   - Final calculated m/z: {sodiated_mass:.5f}")
    
    print("\n---")
    print("Conclusion: All three linkage isomers (a2-3/a2-3, a2-3/a2-6, a2-6/a2-6) have the same elemental composition and thus the same mass.")
    
    return sodiated_mass

# Execute the calculation and store the final answer
final_mass = calculate_glycan_mass()

# The final answer format
print(f"\nFinal Answer: The expected m/z for the singly-sodiated ion for all three glycans is {final_mass:.4f}.")
print(f"<<<{final_mass:.4f}>>>")