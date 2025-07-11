def solve_h2_decomposition():
    """
    This script determines the maximum number of symmetry-adapted Hilbert spaces
    for the H2 molecule in a minimal basis by analyzing its electronic configurations
    and their associated symmetries.
    """
    print("Analyzing the Fock space of the H2 molecule in a minimal basis.")
    print("The minimal basis gives two molecular orbitals: sigma_g (g) and sigma_u (u).")
    print("We consider all possible ways to place 2 electrons into these orbitals.")
    print("-" * 70)

    # We will identify the unique symmetry labels.
    # The full symmetry label is written as ^{2S+1}Λ_{g/u}^{+/-}
    # In our minimal basis, all states are Σ (Lambda=0) and + states.
    # We will use a set to store the unique symmetry labels found.
    symmetries = set()
    
    # --- Configuration 1: (σ_g)² ---
    # Spatial symmetry: g ⊗ g = g
    # Spin: Pauli exclusion principle forces a singlet state (S=0, Multiplicity=1)
    sym1 = "¹Σ_g⁺"
    symmetries.add(sym1)
    print("Configuration (σ_g)²:")
    print(f"  - Spatial Symmetry: g ⊗ g = g")
    print(f"  - Spin Symmetry: Singlet (S=0) due to Pauli exclusion")
    print(f"  - Resulting State Symmetry: {sym1}")
    print("")

    # --- Configuration 2: (σ_u)² ---
    # Spatial symmetry: u ⊗ u = g
    # Spin: Pauli exclusion forces a singlet state (S=0, Multiplicity=1)
    sym2 = "¹Σ_g⁺"
    symmetries.add(sym2) # This will not add a new element to the set
    print("Configuration (σ_u)²:")
    print(f"  - Spatial Symmetry: u ⊗ u = g")
    print(f"  - Spin Symmetry: Singlet (S=0) due to Pauli exclusion")
    print(f"  - Resulting State Symmetry: {sym2}")
    print(f"  - Note: This state belongs to the same Hilbert space as the one from (σ_g)².")
    print("")

    # --- Configuration 3: (σ_g)¹(σ_u)¹ ---
    # Spatial symmetry: g ⊗ u = u
    # Spin: Electrons are in different spatial orbitals, so both singlet and triplet are possible.
    # Singlet state (S=0, Multiplicity=1)
    sym3_singlet = "¹Σ_u⁺"
    symmetries.add(sym3_singlet)
    # Triplet state (S=1, Multiplicity=3)
    sym3_triplet = "³Σ_u⁺"
    symmetries.add(sym3_triplet)
    print("Configuration (σ_g)¹(σ_u)¹:")
    print(f"  - Spatial Symmetry: g ⊗ u = u")
    print(f"  - This configuration gives rise to two distinct symmetry states:")
    print(f"    1. Spin Singlet (S=0) -> State Symmetry: {sym3_singlet}")
    print(f"    2. Spin Triplet (S=1) -> State Symmetry: {sym3_triplet}")
    print("")
    
    print("-" * 70)
    print("Summary of unique symmetry-adapted Hilbert spaces found:")
    # Sort the list for consistent ordering
    sorted_symmetries = sorted(list(symmetries))
    for i, space in enumerate(sorted_symmetries):
        print(f"  {i+1}. {space}")
        
    print("\nThe Fock space F decomposes into the direct sum of these spaces:")
    decomposition_str = " ⊕ ".join([f"H({s})" for s in sorted_symmetries])
    print(f"  F = {decomposition_str}")
    
    print("\nThe maximum number of symmetry-adapted Hilbert spaces is the count of these unique symmetries.")
    
    # Define the number of spaces for each unique symmetry type
    num_Sg1_space = 1 # The ¹Σ_g⁺ space
    num_Su1_space = 1 # The ¹Σ_u⁺ space
    num_Su3_space = 1 # The ³Σ_u⁺ space
    total_spaces = num_Sg1_space + num_Su1_space + num_Su3_space
    
    print(f"\nFinal Calculation: {num_Sg1_space} (for ¹Σ_g⁺) + {num_Su1_space} (for ¹Σ_u⁺) + {num_Su3_space} (for ³Σ_u⁺) = {total_spaces}")
    print(f"\nTherefore, the maximum number of spaces is {total_spaces}.")

solve_h2_decomposition()
<<<3>>>