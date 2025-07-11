def solve_h2_fock_space():
    """
    This function determines the number of symmetry-adapted Hilbert spaces
    for the H2 molecule in a minimal basis by enumerating the unique
    symmetry blocks in its Fock space.
    """
    
    # Each tuple represents a unique symmetry-adapted Hilbert space.
    # The format is (Number of Electrons N, Term Symbol, Dimension of the space)
    # The Term Symbol encodes total spin (superscript), axial angular momentum (Sigma),
    # and inversion symmetry (subscript g/u). All states here are Sigma+.
    subspaces = [
        # N=0 sector (1 state)
        (0, "¹Σg⁺", 1),
        
        # N=1 sector (4 states)
        (1, "²Σg⁺", 2),
        (1, "²Σu⁺", 2),

        # N=2 sector (6 states)
        (2, "¹Σg⁺", 2),
        (2, "³Σu⁺", 3),
        (2, "¹Σu⁺", 1),
        
        # N=3 sector (4 states)
        (3, "²Σg⁺", 2),
        (3, "²Σu⁺", 2),
        
        # N=4 sector (1 state)
        (4, "¹Σg⁺", 1)
    ]

    print("The Fock space of H2 in a minimal basis decomposes into the following symmetry-adapted Hilbert spaces:")
    print("-" * 80)
    print(f"{'Space #':<10} {'N':<5} {'Term Symbol':<15} {'Basis States':<20}")
    print("-" * 80)

    total_states_check = 0
    for i, (n, term, dim) in enumerate(subspaces):
        print(f"{i+1:<10} {n:<5} {term:<15} (Dimension: {dim})")
        total_states_check += dim
        
    print("-" * 80)
    print(f"Check: The sum of dimensions is {total_states_check}, which correctly equals the total number of states (2⁴=16).")
    
    max_subspaces = len(subspaces)
    
    # Generating the "final equation" text
    counts_per_N = {}
    for n, _, _ in subspaces:
        counts_per_N[n] = counts_per_N.get(n, 0) + 1
        
    equation_parts = [str(counts_per_N.get(i, 0)) for i in range(5)]
    equation_str = " + ".join(equation_parts)
    
    print(f"\nNumber of subspaces for N=0, 1, 2, 3, 4 are {equation_parts[0]}, {equation_parts[1]}, {equation_parts[2]}, {equation_parts[3]}, and {equation_parts[4]} respectively.")
    print(f"The total number of symmetry-adapted Hilbert spaces is the sum: {equation_str} = {max_subspaces}")
    print("\nTherefore, the maximum number of symmetry-adapted Hilbert spaces is:")
    print(max_subspaces)

solve_h2_fock_space()
<<<9>>>