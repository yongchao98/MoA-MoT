def solve():
    """
    This script identifies compound C and calculates its molecular weight.
    """
    # Step 1: Explain the deduced structure of compound C
    compound_name = "2,7,12-trihydroxy-4,9-bis(diethylamino)dibenzo[c,h]xanthen-6a-ylium"
    molecular_formula = "C29H30N2O4+"
    
    print(f"The final product, Compound C, is identified as: {compound_name}")
    print(f"Its molecular formula is: {molecular_formula}")
    print("\nCalculating the molecular weight of the cation:")

    # Step 2: Define atomic weights
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Step 3: Define atom counts from the formula C29H30N2O4
    atom_counts = {
        'C': 29,
        'H': 30,
        'N': 2,
        'O': 4
    }

    # Step 4: Calculate and print the molecular weight step-by-step
    c_mass = atom_counts['C'] * atomic_weights['C']
    h_mass = atom_counts['H'] * atomic_weights['H']
    n_mass = atom_counts['N'] * atomic_weights['N']
    o_mass = atom_counts['O'] * atomic_weights['O']
    
    total_mass = c_mass + h_mass + n_mass + o_mass

    print(f"Contribution from Carbon: {atom_counts['C']} * {atomic_weights['C']} = {c_mass:.3f}")
    print(f"Contribution from Hydrogen: {atom_counts['H']} * {atomic_weights['H']} = {h_mass:.3f}")
    print(f"Contribution from Nitrogen: {atom_counts['N']} * {atomic_weights['N']} = {n_mass:.3f}")
    print(f"Contribution from Oxygen: {atom_counts['O']} * {atomic_weights['O']} = {o_mass:.3f}")
    
    print("-" * 30)
    print(f"Total Molecular Weight = {c_mass:.3f} + {h_mass:.3f} + {n_mass:.3f} + {o_mass:.3f} = {total_mass:.3f} g/mol")

solve()