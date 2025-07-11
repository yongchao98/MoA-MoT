def calculate_modification_mass():
    """
    Calculates the mass of the variable modification 'x' on cysteine.

    The plan is as follows:
    1. Identify the chemical probe: 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid.
    2. Determine its reaction with cysteine (Michael addition).
    3. Model the cleavage by formic acid, which breaks the amide bond in the probe.
    4. Identify the fragment that remains attached to the cysteine.
    5. Calculate the mass of this remaining fragment.
    """

    # Atomic masses (using integer values for clarity, as per the question's context)
    atomic_mass = {'C': 12, 'H': 1, 'O': 16, 'N': 14}

    # Step 1: Define the initial probe that reacts with Cysteine.
    # The structure is HOOC-C(=CH2)-CH2-C(=O)-NH-CH2-C≡CH
    # Formula: C9H9NO4
    probe_formula = {'C': 9, 'H': 9, 'N': 1, 'O': 4}
    probe_mass = sum(atomic_mass[atom] * count for atom, count in probe_formula.items())
    
    # Step 2 & 3: After Michael addition, the probe is cleaved by formic acid at the amide bond.
    # The bond ...-C(=O)-|-NH-... is broken.

    # Step 4: Identify the fragment left on cysteine.
    # The structure of the remaining fragment is -CH2-CH(COOH)-CH2-C(=O)-
    # Let's determine its formula.
    # C: 1(CH2)+1(CH)+1(COOH)+1(CH2)+1(C=O) = 5
    # H: 2(CH2)+1(CH)+1(COOH)+2(CH2) = 6
    # O: 2(COOH)+1(C=O) = 3
    final_adduct_formula = {'C': 5, 'H': 6, 'O': 3}

    # Step 5: Calculate the mass of the final adduct. This is the value 'x'.
    x_mass = sum(atomic_mass[atom] * count for atom, count in final_adduct_formula.items())
    
    # Let's also identify the lost fragment to check our work.
    # The lost fragment is -NH-CH2-C≡CH. Its formula is C3H4N.
    # However, the cleavage is a fragmentation, not hydrolysis, based on our derived mass.
    # Let's show the calculation clearly.
    
    c_count = final_adduct_formula['C']
    h_count = final_adduct_formula['H']
    o_count = final_adduct_formula['O']
    
    c_mass = c_count * atomic_mass['C']
    h_mass = h_count * atomic_mass['H']
    o_mass = o_count * atomic_mass['O']
    
    print("The final modification 'x' is the mass of the chemical group that remains on the cysteine after cleavage.")
    print("The chemical formula of this group is C{}H{}O{}.".format(c_count, h_count, o_count))
    print("\nCalculating the mass of this group:")
    print("Mass = (Number of Carbons * Mass of C) + (Number of Hydrogens * Mass of H) + (Number of Oxygens * Mass of O)")
    print("Mass = ({} * {}) + ({} * {}) + ({} * {})".format(c_count, atomic_mass['C'], h_count, atomic_mass['H'], o_count, atomic_mass['O']))
    print("Mass = {} + {} + {}".format(c_mass, h_mass, o_mass))
    print("Mass of x = {}".format(x_mass))

calculate_modification_mass()