import collections

def get_compound_c_properties():
    """
    Calculates and prints the molecular properties of Compound C.
    """
    # Structure of C: tris(2,4,6-trihydroxyphenyl)(diethylamino)methane
    # Chemical name: 4,4',4''-((diethylamino)methanetriyl)tris(benzene-1,3,5-triol)
    
    # 1. Atom counting
    # Central carbon
    central_carbon = {'C': 1}
    # Diethylamino group: -N(CH2CH3)2
    diethylamino_group = {'C': 4, 'H': 10, 'N': 1}
    # One 2,4,6-trihydroxyphenyl group: -C6H2(OH)3
    # Ring: C6, H2. Three OH groups: O3, H3. Total: C6, H5, O3
    # The bond to the central carbon replaces one H on the ring. So, one group is C6H2(OH)3
    # Let's count atoms per group: Ring carbons = 6, Ring hydrogens = 2. Hydroxyls = 3*(O+H) = 3O, 3H.
    # Total per phenyl group: C=6, H=5, O=3
    # Wait, the H count is C6H2(OH)3 -> C6H5O3.
    # Let's verify: Benzene is C6H6. Phenyl is C6H5. Attaching 3 OH groups gives C6H2(OH)3.
    # Ring hydrogens: 6 (benzene) - 1 (attachment) - 3 (for OH) = 2. So C6H2.
    # Hydroxyl hydrogens: 3. Total H = 2 + 3 = 5. So one group is C6H5O3.
    trihydroxyphenyl_group = {'C': 6, 'H': 5, 'O': 3}
    
    # Total atom count
    # Three trihydroxyphenyl groups
    num_phenyl_groups = 3
    
    total_atoms = collections.Counter()
    total_atoms.update(central_carbon)
    total_atoms.update(diethylamino_group)
    total_atoms.update(collections.Counter(trihydroxyphenyl_group) * num_phenyl_groups)
    
    # 2. Atomic weights (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }
    
    # 3. Molecular Formula and Weight Calculation
    c_atoms = total_atoms['C']
    h_atoms = total_atoms['H']
    n_atoms = total_atoms['N']
    o_atoms = total_atoms['O']
    
    molecular_formula = f"C{c_atoms}H{h_atoms}N{n_atoms}O{o_atoms}"
    
    mw = (c_atoms * atomic_weights['C'] + 
          h_atoms * atomic_weights['H'] + 
          n_atoms * atomic_weights['N'] + 
          o_atoms * atomic_weights['O'])

    # 4. Print results
    print(f"The final product is Compound C: tris(2,4,6-trihydroxyphenyl)(diethylamino)methane.")
    print("-" * 30)
    print("Determining Molecular Formula and Weight:")
    print(f"Number of Carbon (C) atoms: {c_atoms}")
    print(f"Number of Hydrogen (H) atoms: {h_atoms}")
    print(f"Number of Nitrogen (N) atoms: {n_atoms}")
    print(f"Number of Oxygen (O) atoms: {o_atoms}")
    print("-" * 30)
    print(f"Molecular Formula: {molecular_formula}")
    print("-" * 30)
    print("Calculation of Molecular Weight (g/mol):")
    # Output each number in the final equation
    print(f"MW = ({c_atoms} * {atomic_weights['C']}) + ({h_atoms} * {atomic_weights['H']}) + ({n_atoms} * {atomic_weights['N']}) + ({o_atoms} * {atomic_weights['O']})")
    print(f"MW = {c_atoms * atomic_weights['C']:.3f} + {h_atoms * atomic_weights['H']:.3f} + {n_atoms * atomic_weights['N']:.3f} + {o_atoms * atomic_weights['O']:.3f}")
    print(f"MW = {mw:.3f} g/mol")

get_compound_c_properties()