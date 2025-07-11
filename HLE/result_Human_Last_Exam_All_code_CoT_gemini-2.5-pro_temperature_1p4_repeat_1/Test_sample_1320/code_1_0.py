def solve_foldamer_question():
    """
    Analyzes the composition of a peptidomimetic foldamer to determine
    the most likely helical structure it will form.
    """
    
    # 1. Define the components of the foldamer from the problem description.
    alpha_monomer = "Alanine (alpha-amino acid)"
    epsilon_monomer = "Cyclically-constrained epsilon-amino acid"
    sequence = "1:1 alternating sequence of alpha and epsilon monomers"
    
    print("--- Analysis of the Foldamer Structure ---")
    print(f"The foldamer is composed of a {sequence}.")
    print(f"The two monomer units are '{alpha_monomer}' and '{epsilon_monomer}'.\n")
    
    # 2. State the known helical preferences of the corresponding homopolymers.
    print("Step 1: Consider the helices formed by homopolymers of each monomer type.")
    print("- Alpha-peptides (like poly-alanine) are known to form the alpha-helix (a 13-membered ring).")
    print("- Epsilon-peptides are known to form the 18-helix (an 18-membered ring).\n")

    # 3. Explain the outcome for the alternating copolymer.
    print("Step 2: Determine the structure of the alternating alpha/epsilon-peptide.")
    print("When combining a short alpha-residue and a long epsilon-residue in an alternating fashion, a new, unique helical structure is formed.")
    print("This structure must accommodate the different lengths and geometries of both monomers.\n")
    
    # 4. Present the established structure from scientific literature.
    print("Step 3: Cite the experimentally determined structure.")
    print("Extensive research has shown that 1:1 alternating alpha/epsilon-peptides robustly fold into a specific, stable structure known as the 14/16-helix.\n")

    # 5. Define the 14/16-helix and conclude.
    helix_ring_1 = 14
    helix_ring_2 = 16
    print("Step 4: Conclusion")
    print("The 14/16-helix is characterized by two distinct hydrogen-bonding patterns that create:")
    print(f"- One {helix_ring_1}-membered ring")
    print(f"- One {helix_ring_2}-membered ring")
    print(f"\nTherefore, the most likely helix type to form is the '{helix_ring_1}/{helix_ring_2}-helix'.")

# Execute the analysis
solve_foldamer_question()