import itertools
from collections import defaultdict

def solve_h2_fock_space():
    """
    Decomposes the Fock space of H2 in a minimal basis into its
    symmetry-adapted Hilbert spaces and counts them.
    """
    # Step 1: Define the basis spin-orbitals
    # (spatial_parity, spin_projection)
    # g = gerade, u = ungerade
    spin_orbitals = [
        ('g', 0.5),   # sigma_g(alpha)
        ('g', -0.5),  # sigma_g(beta)
        ('u', 0.5),   # sigma_u(alpha)
        ('u', -0.5),  # sigma_u(beta)
    ]
    num_spin_orbitals = len(spin_orbitals)

    # This set will store the unique symmetry blocks found
    # Each block is defined by a tuple (N, S, Parity)
    symmetry_blocks = set()

    print("Analyzing the H2 molecule in a minimal basis.")
    print("The basis consists of 4 spin-orbitals: sigma_g(alpha), sigma_g(beta), sigma_u(alpha), sigma_u(beta).")
    print("The Fock space is decomposed into symmetry-adapted Hilbert spaces, labeled by (N, S, Parity):")
    print("  N = Number of electrons")
    print("  S = Total spin quantum number")
    print("  Parity = Spatial inversion symmetry ('g' for gerade, 'u' for ungerade)\n")
    print("The direct sum of these spaces forms the Fock space. Let's find them:")
    print("-" * 60)

    # Step 2: Iterate through all possible numbers of electrons (sectors of Fock space)
    for N in range(num_spin_orbitals + 1):
        # This nested dictionary will store the counts of states for each M_S value,
        # grouped by their spatial parity.
        # Format: {parity: {M_S: count}}
        states_by_parity_ms = defaultdict(lambda: defaultdict(int))

        # Step 3: Generate all possible configurations (Slater determinants) for N electrons
        for config in itertools.combinations(spin_orbitals, N):
            # Calculate M_S (total spin projection)
            ms = sum(so[1] for so in config)

            # Calculate spatial parity
            num_u_orbitals = sum(1 for so in config if so[0] == 'u')
            parity = 'g' if num_u_orbitals % 2 == 0 else 'u'

            states_by_parity_ms[parity][ms] += 1

        # Step 4: For each (N, Parity) group, deduce the S multiplets
        for parity, ms_counts in states_by_parity_ms.items():
            # This is the core of the spin-decomposition algorithm.
            # The number of multiplets of spin S, d(S), can be found from the
            # number of states with M_S=S, c(S), via the formula: d(S) = c(S) - c(S+1)
            
            # Get a sorted list of unique non-negative spin values present
            possible_S_values = sorted([s for s in ms_counts.keys() if s >= 0], reverse=True)

            for S in possible_S_values:
                # Get the number of states with M_S = S and M_S = S+1
                count_S = ms_counts.get(S, 0)
                count_S_plus_1 = ms_counts.get(S + 1, 0)

                # Calculate the number of multiplets with total spin S
                num_multiplets_S = count_S - count_S_plus_1

                if num_multiplets_S > 0:
                    # We found a new symmetry block (or multiple identical ones).
                    # A block is uniquely defined by (N, S, Parity).
                    block = (N, S, parity)
                    symmetry_blocks.add(block)
                    
                    # We need to subtract the contribution of the found multiplet(s)
                    # from the ms_counts to not double-count them for lower S values.
                    for i in range(num_multiplets_S):
                        for ms_val in [s for s in [S, S-1, S-2, -S] if abs(s) <= S]:
                            if S-int(S) == 0.5: # Handle half-integer spins
                                ms_val_range = [S - j for j in range(int(2*S+1))]
                            else: # Handle integer spins
                                ms_val_range = range(-int(S), int(S) + 1)

                        for ms_val in ms_val_range:
                           ms_counts[ms_val] -= 1


    # Step 5: Print the results
    print("Found the following symmetry-adapted Hilbert spaces (N, S, Parity):")
    sorted_blocks = sorted(list(symmetry_blocks))
    for block in sorted_blocks:
        print(f"  {block}")
    
    print("-" * 60)
    print(f"The equation for the Fock space decomposition is a direct sum of these {len(sorted_blocks)} spaces:")
    print("F = " + " ⊕ ".join([f"H({b[0]},{b[1]},'{b[2]}')" for b in sorted_blocks]))
    print("-" * 60)
    
    max_num_spaces = len(symmetry_blocks)
    print(f"The maximum number of symmetry-adapted Hilbert spaces that the Fock space of H2 can be decomposed into is: {max_num_spaces}")
    
    return max_num_spaces

if __name__ == '__main__':
    solve_h2_fock_space()