def rank_epitopes():
    """
    Analyzes and ranks peptide epitopes based on their predicted binding
    affinity to the H2-Kd MHC molecule, explaining the reasoning step-by-step.
    """
    # The epitope sequences for ranking.
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    print("Step 1: Understand the binding rules for the H2-Kd MHC allele.")
    print("H2-Kd binds 9-amino acid peptides. Binding affinity is determined by key anchor residues:")
    print("  - Position 2 (P2): Prefers aromatic residues. Tyrosine (Y) is optimal, and Phenylalanine (F) is also very good.")
    print("  - Position 9 (P9): Must be a hydrophobic residue. Leucine (L), Isoleucine (I), and Valine (V) are good anchors.\n")

    print("Step 2: Evaluate each epitope based on these rules, using E1 as the high-affinity baseline.\n")

    # --- Analysis of each Epitope ---
    print("--- Analysis of E3 (TFQRTRALK) ---")
    print("P2 Anchor is F (good).")
    print("P9 Anchor is K (Lysine). Lysine is positively charged and not hydrophobic.")
    print("Conclusion: A non-permissive P9 anchor drastically reduces binding. E3 will be the WEAKEST binder.\n")

    print("--- Analysis of E2 (TFQRTRALV) vs. E1 (TYQRTRALV) ---")
    print("E1 has the optimal P2 anchor (Y). E2 has a slightly less optimal, but still strong, P2 anchor (F).")
    print("Both have the same good P9 anchor (V) and the same core sequence.")
    print("Conclusion: The optimal P2 anchor makes E1 a stronger binder than E2.\n")

    print("--- Analysis of E4 (TYQRMFALV) vs. E1 (TYQRTRALV) ---")
    print("Both have identical, optimal P2 (Y) and P9 (V) anchors.")
    print("E4's sequence differs from E1 at P5 (T -> M) and P6 (R -> F).")
    print("These changes increase the hydrophobicity of the peptide core, which generally enhances interaction with the MHC binding groove.")
    print("Conclusion: E4 is predicted to be a STRONGER binder than E1.\n")
    
    print("--- Analysis of E5 (TYIPMRALV) vs. E1 (TYQRTRALV) ---")
    print("Both have identical, optimal P2 (Y) and P9 (V) anchors.")
    print("E5 differs from E1 at P3 (Q -> I), P4 (R -> P), and P5 (T -> M).")
    print("  - P3 (Q -> I) and P5 (T -> M) are favorable changes to better secondary anchors.")
    print("  - P4 (R -> P) introduces a Proline, which is structurally disruptive and generally unfavorable for binding.")
    print("Conclusion: The strong positive effects at P3 and P5 are thought to outweigh the negative effect of the P4 Proline, leading to a net binding affinity greater than E1, but the disruption makes it weaker than E4.\n")

    print("Step 3: Assemble the final ranking from highest to lowest affinity.")
    # The final ranking is synthesized from the pairwise comparisons.
    final_ranking = ["E4", "E5", "E1", "E2", "E3"]
    print("1st (Highest Affinity): E4, due to optimal anchors and improved core residues.")
    print("2nd: E5, due to strong positive changes outweighing a negative one.")
    print("3rd: E1, the original high-affinity baseline.")
    print("4th: E2, strong but with a slightly suboptimal P2 anchor compared to E1.")
    print("5th (Lowest Affinity): E3, due to the detrimental P9 anchor.\n")

    print("Final Rank Order: E4, E5, E1, E2, E3")

# Execute the analysis and print the result.
rank_epitopes()
<<<B>>>