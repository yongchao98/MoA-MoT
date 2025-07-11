def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to the H-2Kd MHC allele.
    """
    # Define the epitopes
    epitopes = {
        "E1": "T Y Q R T R A L V",
        "E2": "T F Q R T R A L V",
        "E3": "T F Q R T R A L K",
        "E4": "T Y Q R M F A L V",
        "E5": "T Y I P M R A L V"
    }

    # The ranking is determined by analyzing the anchor residues (P2 and P9)
    # and other mutations based on known H-2Kd binding motifs.
    # H-2Kd anchors: P2=Y/F, P9=V/L/I.
    #
    # 1. E1 (TYQRTRALV): P2=Y, P9=V. The high-affinity reference.
    # 2. E4 (TYQRMFALV): P2=Y, P9=V. Optimal anchors, non-disruptive changes elsewhere. High affinity.
    # 3. E5 (TYIPMRALV): P2=Y, P9=V. Optimal anchors, but a disruptive Proline (P) at P4 reduces affinity.
    # 4. E2 (TFQRTRALV): P2=F (good), P9=V. A change at an anchor position from the reference makes it weaker than the top 3.
    # 5. E3 (TFQRTRALK): P2=F (good), P9=K (poor). A non-conservative change at the critical P9 anchor drastically reduces affinity.
    
    ranked_order = ["E1", "E4", "E5", "E2", "E3"]
    
    print("Ranking of epitopes from highest to lowest expected amount complexed with H2-Kd:")
    for i, epi_id in enumerate(ranked_order):
        print(f"{i+1}. {epi_id}: {epitopes[epi_id]}")

# Execute the function to print the ranked list.
rank_epitopes()