import collections

def solve_h2_fock_space():
    """
    Analyzes the Fock space of the H2 molecule in a minimal basis to find the
    maximum number of symmetry-adapted Hilbert spaces.
    """
    # Step 1: Define the 6 Configuration State Functions (CSFs)
    # These form a basis for the 6-dimensional Fock space. Each is defined by its
    # quantum numbers: S (total spin), M_S (spin projection), and Parity.
    csfs = [
        # CSF 1: From the (sigma_g)^2 configuration. This is the main component of the ground state.
        {'name': '¹Σ_g⁺(1)', 'S': 0, 'M_S': 0, 'Parity': 'g'},
        # CSF 2: From the (sigma_u)^2 configuration. This is a doubly-excited state.
        {'name': '¹Σ_g⁺(2)', 'S': 0, 'M_S': 0, 'Parity': 'g'},
        # CSF 3: The singlet state from the (sigma_g)^1(sigma_u)^1 configuration. A singly-excited state.
        {'name': '¹Σ_u⁺', 'S': 0, 'M_S': 0, 'Parity': 'u'},
        # CSFs 4, 5, 6: The three components of the triplet state from the (sigma_g)^1(sigma_u)^1 configuration.
        {'name': '³Σ_u⁺ (M_S=+1)', 'S': 1, 'M_S': 1, 'Parity': 'u'},
        {'name': '³Σ_u⁺ (M_S=0)', 'S': 1, 'M_S': 0, 'Parity': 'u'},
        {'name': '³Σ_u⁺ (M_S=-1)', 'S': 1, 'M_S': -1, 'Parity': 'u'},
    ]

    # Step 2: Decompose the Fock space using symmetries
    # We group CSFs by their (S, M_S, Parity) quantum numbers. The Hamiltonian is
    # block-diagonal with respect to this grouping. Each group forms a basis for a
    # symmetry-adapted Hilbert space.
    symmetry_spaces = collections.defaultdict(list)
    for csf in csfs:
        # Create a unique key for each set of quantum numbers
        key = (f"S={csf['S']}", f"M_S={csf['M_S']}", f"Parity='{csf['Parity']}'")
        symmetry_spaces[key].append(csf['name'])

    # Step 3: Display the results
    print("The 6-dimensional Fock space of H2 in a minimal basis decomposes into the following symmetry-adapted Hilbert spaces:")
    print("-" * 80)

    space_dimensions = []
    for i, (key, states) in enumerate(symmetry_spaces.items()):
        dimension = len(states)
        space_dimensions.append(str(dimension))
        print(f"Space {i+1}:")
        print(f"  Quantum Numbers: {', '.join(key)}")
        print(f"  Dimension: {dimension}")
        print(f"  Basis States: {', '.join(states)}")
        print("-" * 80)

    # The total number of spaces is the number of unique groups we found.
    max_spaces = len(symmetry_spaces)

    print("The Fock space is a direct sum of these subspaces.")
    print("The total dimension is 6, which is the sum of the dimensions of the individual subspaces:")
    # We print the decomposition equation as requested
    print(f"6 = {' + '.join(sorted(space_dimensions, reverse=True))}")
    print("\nThis decomposition fully block-diagonalizes the Hamiltonian matrix.")
    print("The maximum number of symmetry-adapted Hilbert spaces the Fock space can be decomposed into is therefore:")
    print(max_spaces)

solve_h2_fock_space()