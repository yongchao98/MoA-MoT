import collections

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to the H2-Kd MHC allele.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }

    print("--- Epitope Binding Analysis for H2-Kd ---")
    print("The binding affinity of a peptide to H2-Kd is primarily determined by anchor residues at Position 2 (P2) and Position 9 (P9).")
    print("Ideal P2 anchor: Y (Tyrosine)")
    print("Ideal P9 anchor: V (Valine), L (Leucine), or I (Isoleucine)\n")

    # Step 1: Analyze the P9 C-terminal anchor. This is the most critical position.
    print("Step 1: Evaluating the C-terminal anchor (P9)...")
    weakest_binder = None
    for name, seq in epitopes.items():
        p9_residue = seq[8]
        if p9_residue not in ['V', 'L', 'I']:
            print(f"Epitope {name} ('{seq}') has an unfavorable anchor '{p9_residue}' at P9.")
            print("A charged residue (K) is highly incompatible with the hydrophobic P9 pocket.")
            print(f"Result: {name} will be the weakest binder.\n")
            weakest_binder = name

    remaining_epitopes = {k: v for k, v in epitopes.items() if k != weakest_binder}

    # Step 2: Analyze the P2 anchor for the remaining epitopes.
    print("Step 2: Evaluating the P2 anchor for the remaining epitopes...")
    p2_suboptimal_binder = None
    for name, seq in remaining_epitopes.items():
        p2_residue = seq[1]
        if p2_residue == 'F':
            print(f"Epitope {name} ('{seq}') has a Phenylalanine (F) at P2.")
            print("While F is an acceptable anchor, it is generally suboptimal compared to Tyrosine (Y).")
            print(f"Result: {name} will bind less tightly than epitopes with Y at P2.\n")
            p2_suboptimal_binder = name

    top_candidates = {k: v for k, v in remaining_epitopes.items() if k != p2_suboptimal_binder}
    
    # Step 3: Rank the top candidates based on central residues.
    print("Step 3: Evaluating central residues for top candidates (E1, E4, E5)...")
    print("All have optimal anchors (Y at P2, V at P9). Their rank depends on other mutations.")
    
    e1_seq = top_candidates['E1']
    e4_seq = top_candidates['E4']
    e5_seq = top_candidates['E5']
    
    print(f"E1 ('{e1_seq}') is the reference high-affinity epitope. It is the benchmark for optimal binding.")
    print(f"E5 ('{e5_seq}') contains a Proline ('P') at P4. Proline's rigid structure disrupts the peptide backbone, likely reducing affinity significantly.")
    print(f"E4 ('{e4_seq}') has substitutions at non-anchor positions (P5, P6). These changes are less disruptive than the Proline in E5.")
    print("Result: The predicted binding order for this group is E1 > E4 > E5.\n")
    
    final_ranking = ['E1', 'E4', 'E5', p2_suboptimal_binder, weakest_binder]

    # Final conclusion
    print("--- Final Ranking ---")
    print("Combining the analyses gives the final rank order from highest to lowest expected amount complexed with H2-Kd:")
    print(", ".join(final_ranking))

rank_epitopes()