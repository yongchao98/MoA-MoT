def solve_h2_fock_space():
    """
    Calculates the number of symmetry-adapted Hilbert spaces for the H2 molecule
    in a minimal basis by enumerating all possible term symbols.
    """
    print("Decomposition of the H2 Molecule's Fock Space in a Minimal Basis")
    print("=================================================================\n")
    print("The electronic Hamiltonian commutes with operators for electron number (N),")
    print("total spin (S), and spatial point group symmetry (Gamma). The Fock space")
    print("can be decomposed into blocks, where each block has a unique set of")
    print("quantum numbers (N, S, Gamma).\n")
    
    # The set will store unique symmetry spaces, identified by a tuple (N, S, Gamma).
    # Gamma is represented as a string like '1Sg+' (for multiplicity 1, Sigma_g+)
    all_spaces = set()
    
    # Store the count of spaces for each N to construct the final sum
    counts_per_n = {}

    # --- N = 0 ---
    n = 0
    # Config: (sigma_g)^0(sigma_u)^0 -> Vacuum state
    # S=0, Parity=g, Symbol: 1Sg+
    spaces_n0 = {(n, 0, 'Sg+')}
    counts_per_n[n] = len(spaces_n0)
    all_spaces.update(spaces_n0)
    print(f"For N={n} (0 electrons):")
    print(f"  - The vacuum state has S=0 and is totally symmetric (1Sg+).")
    print(f"  - Number of symmetry spaces: {counts_per_n[n]}\n")

    # --- N = 1 ---
    n = 1
    # Config: (sigma_g)^1 -> S=1/2, Parity=g -> 2Sg+
    # Config: (sigma_u)^1 -> S=1/2, Parity=u -> 2Su+
    spaces_n1 = {(n, 0.5, 'Sg+'), (n, 0.5, 'Su+')}
    counts_per_n[n] = len(spaces_n1)
    all_spaces.update(spaces_n1)
    print(f"For N={n} (1 electron):")
    print(f"  - Config (sg)^1 gives one space: 2Sg+")
    print(f"  - Config (su)^1 gives another space: 2Su+")
    print(f"  - Number of symmetry spaces: {counts_per_n[n]}\n")

    # --- N = 2 ---
    n = 2
    spaces_n2 = set()
    # Config: (sigma_g)^2 -> S=0, Parity=g -> 1Sg+
    spaces_n2.add((n, 0, 'Sg+'))
    # Config: (sigma_u)^2 -> S=0, Parity=g -> 1Sg+ (same symmetry)
    spaces_n2.add((n, 0, 'Sg+'))
    # Config: (sigma_g)^1(sigma_u)^1 -> gives two states
    # S=0, Parity=u -> 1Su+
    spaces_n2.add((n, 0, 'Su+'))
    # S=1, Parity=u -> 3Su+
    spaces_n2.add((n, 1, 'Su+'))
    counts_per_n[n] = len(spaces_n2)
    all_spaces.update(spaces_n2)
    print(f"For N={n} (2 electrons):")
    print(f"  - Configs (sg)^2 and (su)^2 both contribute to the 1Sg+ space.")
    print(f"  - Config (sg)^1(su)^1 gives rise to two distinct spaces: 1Su+ and 3Su+.")
    print(f"  - Total unique spaces are: { {s[2] for s in spaces_n2} }.")
    print(f"  - Number of symmetry spaces: {counts_per_n[n]}\n")

    # --- N = 3 ---
    n = 3
    # Particle-hole symmetry with N=1
    # Config: (sigma_g)^2(sigma_u)^1 -> S=1/2, Parity=u -> 2Su+
    # Config: (sigma_g)^1(sigma_u)^2 -> S=1/2, Parity=g -> 2Sg+
    spaces_n3 = {(n, 0.5, 'Sg+'), (n, 0.5, 'Su+')}
    counts_per_n[n] = len(spaces_n3)
    all_spaces.update(spaces_n3)
    print(f"For N={n} (3 electrons, hole picture of N=1):")
    print(f"  - A hole in the sg orbital gives a 2Sg+ state.")
    print(f"  - A hole in the su orbital gives a 2Su+ state.")
    print(f"  - Number of symmetry spaces: {counts_per_n[n]}\n")

    # --- N = 4 ---
    n = 4
    # Config: (sigma_g)^2(sigma_u)^2 -> Fully occupied
    # S=0, Parity=g -> 1Sg+
    spaces_n4 = {(n, 0, 'Sg+')}
    counts_per_n[n] = len(spaces_n4)
    all_spaces.update(spaces_n4)
    print(f"For N={n} (4 electrons):")
    print(f"  - The filled-shell state is a single 1Sg+ space.")
    print(f"  - Number of symmetry spaces: {counts_per_n[n]}\n")

    # --- Final Calculation ---
    total_count = len(all_spaces)
    sum_parts = [str(counts_per_n[i]) for i in sorted(counts_per_n.keys())]
    sum_expression = " + ".join(sum_parts)
    
    print("-----------------------------------------------------------------")
    print("Final Result:")
    print("The total number of symmetry-adapted Hilbert spaces is the sum of")
    print("the counts from each N-electron sector:")
    print(f"\nTotal Spaces = (N=0) + (N=1) + (N=2) + (N=3) + (N=4)")
    print(f"             = {sum_expression}")
    print(f"             = {total_count}")
    print(f"\nThe maximum number of symmetry-adapted Hilbert spaces is {total_count}.")


if __name__ == "__main__":
    solve_h2_fock_space()
    print("\n<<<9>>>")