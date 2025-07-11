import collections

def calculate_mass(formula, atomic_masses):
    """Calculates the monoisotopic mass of a molecule from its formula."""
    mass = 0.0
    for atom, count in formula.items():
        mass += count * atomic_masses[atom]
    return mass

def main():
    """
    Calculates the m/z of a derivatized glycan.
    """
    # Monoisotopic atomic masses
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Step 1: Base glycan (A2G2S2, a.k.a. Neu5Ac2Gal2GlcNAc4Man3) formula
    # Composition: 4x GlcNAc, 2x Gal, 3x Man, 2x Neu5Ac
    # Formula for a released N-glycan with one reducing end: C84 H138 N6 O62
    base_formula = collections.defaultdict(int, {'C': 84, 'H': 138, 'N': 6, 'O': 62})
    base_mass = calculate_mass(base_formula, atomic_masses)
    
    print("Step 1: Start with the native A2G2S2 glycan.")
    print(f"  - Initial Formula: C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}")
    print(f"  - Initial Mass: {base_mass:.5f} Da\n")

    # Step 2: Amidation
    # Converts 2 COOH groups to CONH2. Change per group: -O, +N, +H
    num_amidations = 2
    amidated_formula = base_formula.copy()
    amidated_formula['O'] -= num_amidations
    amidated_formula['N'] += num_amidations
    amidated_formula['H'] += num_amidations
    amidated_mass = calculate_mass(amidated_formula, atomic_masses)

    print("Step 2: Perform amidation on the two sialic acids (2x COOH -> CONH2).")
    print(f"  - Formula after Amidation: C{amidated_formula['C']}H{amidated_formula['H']}N{amidated_formula['N']}O{amidated_formula['O']}")
    print(f"  - Mass after Amidation: {amidated_mass:.5f} Da\n")

    # Step 3: Permethylation
    # Methylates all free -OH and N-acetyl -NH groups.
    # Total sites = 28 (-OH) + 6 (-NH) = 34 sites.
    # Change per site: replace H with CH3. Net change: +C, +2H.
    num_methylations = 34
    final_formula = amidated_formula.copy()
    final_formula['C'] += num_methylations
    # H change: -num_methylations (from OH/NH) + 3*num_methylations (from CH3) = +2*num_methylations
    final_formula['H'] += 2 * num_methylations
    final_mass = calculate_mass(final_formula, atomic_masses)
    
    print("Step 3: Perform permethylation on all 34 acidic protons (28 -OH and 6 N-acetyl -NH).")
    print(f"  - Final Derivatized Formula: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")
    print(f"  - Final Mass (M): {final_mass:.5f} Da\n")

    # Step 4: Add Sodium for [M+Na]+ ion
    sodiated_mass = final_mass + atomic_masses['Na']

    print("Step 4: Calculate the m/z for the singly sodiated ion [M+Na]⁺.")
    print(f"  - Final m/z: {final_mass:.5f} + {atomic_masses['Na']:.5f} = {sodiated_mass:.5f}\n")

    print("The three glycans are linkage isomers with identical chemical formulas and numbers of reactive sites.")
    print("Therefore, they will all have the same mass after derivatization.")

    print("\n--- Final Results ---")
    print(f"The expected m/z for A2G(4)2S(3)2 is: {sodiated_mass:.4f}")
    print(f"The expected m/z for A2G(4)S(3)S(6) is: {sodiated_mass:.4f}")
    print(f"The expected m/z for A2G(4)2S(6)2 is: {sodiated_mass:.4f}\n")
    
    print("Final mass calculation breakdown:")
    c, h, n, o = final_formula['C'], final_formula['H'], final_formula['N'], final_formula['O']
    c_mass, h_mass, n_mass, o_mass, na_mass = atomic_masses['C'], atomic_masses['H'], atomic_masses['N'], atomic_masses['O'], atomic_masses['Na']
    print(f"({c} * {c_mass:.6f}) + ({h} * {h_mass:.6f}) + ({n} * {n_mass:.6f}) + ({o} * {o_mass:.6f}) + {na_mass:.6f} = {sodiated_mass:.4f}")

if __name__ == '__main__':
    main()
<<<2720.3369>>>