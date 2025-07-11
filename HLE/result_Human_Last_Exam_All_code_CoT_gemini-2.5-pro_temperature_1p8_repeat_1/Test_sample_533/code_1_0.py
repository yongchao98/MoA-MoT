import collections

def analyze_h2_fock_space():
    """
    Analyzes the Fock space of H2 in a minimal basis to find the number of
    symmetry-adapted Hilbert spaces.
    """
    print("### Analysis of H2 Fock Space Decomposition ###\n")
    print("Step 1: Define the basis and symmetries.")
    print("Minimal basis -> Molecular Orbitals: σ_g (gerade), σ_u (ungerade).")
    print("Symmetries used for decomposition:")
    print("  - N: Total number of electrons")
    print("  - S: Total spin quantum number")
    print("  - Γ: Spatial symmetry (irreducible representation of the point group)\n")

    # A list to store the unique Hilbert spaces found, defined by (N, S, Γ)
    # The term symbol '^(2S+1)Γ' is used for the spatial symmetry label Γ.
    hilbert_spaces = []

    # Use a set to track unique spaces to avoid duplicates, e.g., the two CSFs
    # of ¹Σ_g⁺ symmetry in the N=2 sector belong to the same Hilbert space.
    unique_space_labels = set()

    # --- N = 0 Sector ---
    print("--- Analyzing N=0 Sector ---")
    n = 0
    s = 0
    multiplicity = 2 * s + 1
    gamma = "Σ_g⁺"
    term_symbol = f"¹{gamma}"
    space_label = (n, s, term_symbol)
    if space_label not in unique_space_labels:
        print(f"Configuration '()' (vacuum) gives a Hilbert space with symmetry: {term_symbol}")
        unique_space_labels.add(space_label)
        hilbert_spaces.append({'N': n, 'Term': term_symbol, 'Description': 'Vacuum state'})
    print("-" * 30)

    # --- N = 1 Sector ---
    print("--- Analyzing N=1 Sector ---")
    n = 1
    s = 0.5
    multiplicity = int(2 * s + 1)
    
    # Config (σ_g)¹
    gamma_g = "Σ_g⁺"
    term_symbol_g = f"²{gamma_g}"
    space_label_g = (n, s, term_symbol_g)
    if space_label_g not in unique_space_labels:
        print(f"Configuration '(σ_g)¹' gives a space with symmetry: {term_symbol_g}")
        unique_space_labels.add(space_label_g)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_g, 'Description': 'Electron in bonding MO'})
        
    # Config (σ_u)¹
    gamma_u = "Σ_u⁺"
    term_symbol_u = f"²{gamma_u}"
    space_label_u = (n, s, term_symbol_u)
    if space_label_u not in unique_space_labels:
        print(f"Configuration '(σ_u)¹' gives a space with symmetry: {term_symbol_u}")
        unique_space_labels.add(space_label_u)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_u, 'Description': 'Electron in anti-bonding MO'})
    print("-" * 30)

    # --- N = 2 Sector ---
    print("--- Analyzing N=2 Sector ---")
    n = 2
    # Configs (σ_g)² and (σ_u)²
    s_singlet = 0
    gamma_gg = "Σ_g⁺" # from g*g=g and u*u=g
    term_symbol_gg = f"¹{gamma_gg}"
    space_label_gg = (n, s_singlet, term_symbol_gg)
    if space_label_gg not in unique_space_labels:
        print(f"Configurations '(σ_g)²' and '(σ_u)²' together span a space with symmetry: {term_symbol_gg}")
        unique_space_labels.add(space_label_gg)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_gg, 'Description': 'Two closed-shell-like configurations'})

    # Config (σ_g)¹(σ_u)¹
    print(f"Configuration '(σ_g)¹(σ_u)¹' gives rise to two spaces:")
    # Singlet state
    gamma_gu_singlet = "Σ_u⁺" # From symmetric spatial part (g*u=u)
    term_symbol_gu_singlet = f"¹{gamma_gu_singlet}"
    space_label_gu_singlet = (n, s_singlet, term_symbol_gu_singlet)
    if space_label_gu_singlet not in unique_space_labels:
        print(f"  - Singlet state with symmetry: {term_symbol_gu_singlet}")
        unique_space_labels.add(space_label_gu_singlet)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_gu_singlet, 'Description': 'Open-shell singlet state'})
        
    # Triplet state
    s_triplet = 1
    gamma_gu_triplet = "Σ_u⁺" # From anti-symmetric spatial part (g*u=u)
    term_symbol_gu_triplet = f"³{gamma_gu_triplet}"
    space_label_gu_triplet = (n, s_triplet, term_symbol_gu_triplet)
    if space_label_gu_triplet not in unique_space_labels:
        print(f"  - Triplet state with symmetry: {term_symbol_gu_triplet}")
        unique_space_labels.add(space_label_gu_triplet)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_gu_triplet, 'Description': 'Open-shell triplet state'})
    print("-" * 30)

    # --- N = 3 Sector ---
    print("--- Analyzing N=3 Sector (by particle-hole symmetry) ---")
    n = 3
    s = 0.5
    # Corresponds to a single hole in the N=4 configuration (σ_g)²(σ_u)²
    # A hole in σ_g -> same symmetry as an electron in σ_g
    term_symbol_g_hole = "²Σ_g⁺"
    space_label_g_hole = (n, s, term_symbol_g_hole)
    if space_label_g_hole not in unique_space_labels:
        print(f"Configuration '(σ_g)¹(σ_u)²' (a hole in σ_g) has symmetry: {term_symbol_g_hole}")
        unique_space_labels.add(space_label_g_hole)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_g_hole, 'Description': 'Hole in bonding MO'})

    # A hole in σ_u -> same symmetry as an electron in σ_u
    term_symbol_u_hole = "²Σ_u⁺"
    space_label_u_hole = (n, s, term_symbol_u_hole)
    if space_label_u_hole not in unique_space_labels:
        print(f"Configuration '(σ_g)²(σ_u)¹' (a hole in σ_u) has symmetry: {term_symbol_u_hole}")
        unique_space_labels.add(space_label_u_hole)
        hilbert_spaces.append({'N': n, 'Term': term_symbol_u_hole, 'Description': 'Hole in anti-bonding MO'})
    print("-" * 30)

    # --- N = 4 Sector ---
    print("--- Analyzing N=4 Sector ---")
    n = 4
    s = 0
    gamma = "Σ_g⁺"
    term_symbol = f"¹{gamma}"
    space_label = (n, s, term_symbol)
    if space_label not in unique_space_labels:
        print(f"Configuration '(σ_g)²(σ_u)²' (filled shell) has symmetry: {term_symbol}")
        unique_space_labels.add(space_label)
        hilbert_spaces.append({'N': n, 'Term': term_symbol, 'Description': 'Fully filled shell state'})
    print("-" * 30)
    
    # --- Final Summary ---
    num_spaces = len(hilbert_spaces)
    print("\n### Final Result ###")
    print(f"The Fock space is decomposed into the direct sum of {num_spaces} symmetry-adapted Hilbert spaces.")
    print("The decomposition is: F = H₁ ⊕ H₂ ⊕ ... ⊕ H₉, where the spaces Hᵢ are:")
    
    # Sort for deterministic output
    sorted_spaces = sorted(hilbert_spaces, key=lambda x: (x['N'], x['Term']))
    
    for i, space in enumerate(sorted_spaces, 1):
        print(f"  {i}. N={space['N']}, Symmetry={space['Term']:>6s} ({space['Description']})")
        
    print(f"\nThe maximum number of symmetry-adapted Hilbert spaces is {num_spaces}.")


if __name__ == "__main__":
    analyze_h2_fock_space()
    # Final answer as requested by the format
    print("\n<<<9>>>")