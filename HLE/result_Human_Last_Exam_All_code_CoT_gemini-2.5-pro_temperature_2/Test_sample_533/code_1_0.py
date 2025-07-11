import itertools

def analyze_h2_fock_space():
    """
    Analyzes the H2 molecule's Fock space in a minimal basis to determine
    the maximum number of symmetry-adapted Hilbert spaces.
    """
    print("Step-by-step analysis of H2 Fock Space Decomposition")
    print("======================================================")
    print()

    # Step 1: Define the basis
    print("Step 1: Define the basis set")
    print("------------------------------")
    print("The system is the H2 molecule in a minimal basis, which gives two molecular orbitals (MOs):")
    print(" - σ_g (bonding, gerade/symmetric, denoted 'g')")
    print(" - σ_u (antibonding, ungerade/antisymmetric, denoted 'u')")
    print()
    print("Including spin (α for up, β for down), we get 4 spin-orbitals:")
    spin_orbitals = [('g', 'α'), ('g', 'β'), ('u', 'α'), ('u', 'β')]
    for i, so in enumerate(spin_orbitals):
        print(f"  {i+1}. |{so[0]}, {so[1]}>")
    print()

    # Step 2: Construct the Fock space
    print("Step 2: Construct the Fock space")
    print("---------------------------------")
    print("For 2 electrons in 4 spin-orbitals, the number of configurations is C(4, 2) = 6.")
    configs = list(itertools.combinations(spin_orbitals, 2))
    print("The 6 basis configurations are:")
    for i, cfg in enumerate(configs):
        print(f"  Ψ{i+1} = a determinant of |{cfg[0][0]},{cfg[0][1]}> and |{cfg[1][0]},{cfg[1][1]}>")
    print()

    # Step 3: Identify Symmetries
    print("Step 3: Identify Hamiltonian symmetries")
    print("---------------------------------------")
    print("The Hamiltonian commutes with total spin (S²) and its projection (S_z), and spatial point group operators.")
    print("We use these symmetries to classify states. Spatial symmetries combine like: g*g=g, u*u=g, g*u=u.")
    print()

    # Step 4: Decompose the Fock space
    print("Step 4: Decompose Fock space into symmetry-adapted spaces")
    print("---------------------------------------------------------")
    print("We classify the configurations to find the non-mixing subspaces (blocks of the Hamiltonian).")
    print()

    print("--- Subspace 1: ¹Ag ---")
    print("The (σ_g)² and (σ_u)² configurations both have M_s=0 and spatial symmetry 'g'.")
    print("  - Configuration |g,α>|g,β>: S=0 (singlet), M_s=0, Spatial=g -> ¹Ag")
    print("  - Configuration |u,α>|u,β>: S=0 (singlet), M_s=0, Spatial=g -> ¹Ag")
    print("=> These two states form a single 2-dimensional Hilbert space of ¹Ag symmetry.")
    print()

    print("--- Subspace 2: ¹B1u ---")
    print("From the (σ_g)¹(σ_u)¹ configurations with M_s=0, a singlet linear combination can be formed.")
    print("  - State: Combination of |g,α>|u,β> and |g,β>|u,α>")
    print("  - S=0 (singlet), M_s=0, Spatial=g*u=u -> ¹Bu (using a common notation)")
    print("=> This state forms a 1-dimensional Hilbert space of ¹Bu symmetry.")
    print()

    print("--- Subspace 3: ³B1u (M_s = +1) ---")
    print("From the (σ_g)¹(σ_u)¹ configurations, we identify the triplet states.")
    print("  - Configuration |g,α>|u,α>:")
    print("  - S=1 (triplet), M_s=+1, Spatial=g*u=u -> ³Bu")
    print("=> This state forms a 1-dimensional Hilbert space for the M_s=+1 component of the triplet.")
    print()

    print("--- Subspace 4: ³B1u (M_s = 0) ---")
    print("  - State: The triplet linear combination of |g,α>|u,β> and |g,β>|u,α>")
    print("  - S=1 (triplet), M_s=0, Spatial=g*u=u -> ³Bu")
    print("=> This state forms a 1-dimensional Hilbert space for the M_s=0 component of the triplet.")
    print()

    print("--- Subspace 5: ³B1u (M_s = -1) ---")
    print("  - Configuration |g,β>|u,β>:")
    print("  - S=1 (triplet), M_s=-1, Spatial=g*u=u -> ³Bu")
    print("=> This state forms a 1-dimensional Hilbert space for the M_s=-1 component of the triplet.")
    print()


    # Step 5: Final Count
    print("Step 5: Conclusion")
    print("--------------------")
    print("The Fock space is a direct sum of these 5 mutually orthogonal subspaces:")
    print("  Space = (¹Ag space) ⊕ (¹Bu space) ⊕ (³Bu, Ms=+1) ⊕ (³Bu, Ms=0) ⊕ (³Bu, Ms=-1)")
    max_spaces = 5
    print(f"The dimensions of these spaces are 2, 1, 1, 1, and 1, respectively (total dimension is 6).")
    print(f"The maximum number of symmetry-adapted Hilbert spaces is the number of distinct symmetry groups found.")
    print()
    print("Final Answer:")
    print(f"The maximum number of symmetry-adapted Hilbert spaces is {max_spaces}.")

analyze_h2_fock_space()
<<<5>>>