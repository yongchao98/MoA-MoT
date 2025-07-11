def calculate_helix_parameters():
    """
    Calculates the theoretical H-bond ring sizes for an alternating
    alpha/epsilon peptide helix based on i->i+4 interactions.
    """
    # Number of backbone atoms for each residue type
    n_b_alpha = 3  # For Alanine: -N-Ca-C'-
    
    # Model 1: Linear epsilon-amino acid
    n_b_epsilon_linear = 7 # For linear e-AA: -N-Ce-Cd-Cg-Cb-Ca-C'-
    
    # Model 2: Cyclically-constrained epsilon-amino acid (hypothesized shorter path)
    n_b_epsilon_constrained = 6 # e.g., similar to a delta-amino acid

    print("Analysis based on C=O(i) -> H-N(i+4) hydrogen bonds:")
    print("-" * 50)

    # --- Case 1: H-bond between two Alanine residues ---
    # Formula for ring size: 6 + 2 * (number of backbone atoms in epsilon residue)
    print("H-bond pattern: Ala(i) -> Ala(i+4)")
    print("This H-bond must cross two intervening epsilon residues.")
    print("Ring Size = 6 + 2 * n_b(e)\n")

    # Using linear epsilon-AA model
    ring_size_A_linear = 6 + 2 * n_b_epsilon_linear
    print("For a LINEAR epsilon-AA (n_b(e) = 7):")
    print(f"Ring Size = 6 + 2 * {n_b_epsilon_linear} = {ring_size_A_linear} (C{ring_size_A_linear} helix)")
    
    # Using constrained epsilon-AA model
    ring_size_A_constrained = 6 + 2 * n_b_epsilon_constrained
    print("\nFor a CONSTRAINED epsilon-AA (n_b(e) = 6):")
    print(f"Ring Size = 6 + 2 * {n_b_epsilon_constrained} = {ring_size_A_constrained} (C{ring_size_A_constrained} helix)")
    print("-" * 50)

    # --- Case 2: H-bond between two Epsilon-amino acid residues ---
    # Formula for ring size: 9 + (number of backbone atoms in epsilon residue)
    print("H-bond pattern: eAA(i) -> eAA(i+4)")
    print("This H-bond must cross one intervening epsilon residue and two alanine residues.")
    print("Ring Size = 9 + n_b(e)\n")
    
    # Using linear epsilon-AA model
    ring_size_E_linear = 9 + n_b_epsilon_linear
    print("For a LINEAR epsilon-AA (n_b(e) = 7):")
    print(f"Ring Size = 9 + {n_b_epsilon_linear} = {ring_size_E_linear} (C{ring_size_E_linear} helix)")
    
    # Using constrained epsilon-AA model
    ring_size_E_constrained = 9 + n_b_epsilon_constrained
    print("\nFor a CONSTRAINED epsilon-AA (n_b(e) = 6):")
    print(f"Ring Size = 9 + {n_b_epsilon_constrained} = {ring_size_E_constrained} (C{ring_size_E_constrained} helix)")
    print("-" * 50)
    
    print("\nSummary of potential i->i+4 ring sizes:")
    print(f"Linear eAA model: C{ring_size_E_linear} and C{ring_size_A_linear}")
    print(f"Constrained eAA model: C{ring_size_E_constrained} and C{ring_size_A_constrained}")

calculate_helix_parameters()