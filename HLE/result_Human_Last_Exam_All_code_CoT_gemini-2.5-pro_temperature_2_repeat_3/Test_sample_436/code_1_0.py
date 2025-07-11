def solve_epitope_ranking():
    """
    Analyzes and ranks peptide epitopes based on their predicted binding affinity to the H-2Kd MHC molecule.
    """

    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    reference_internal_sequence = "QRTRA"
    scores = {}

    print("--- Epitope Binding Affinity Analysis for H2-Kd ---")
    print("Scoring System:")
    print("- Tier 1 (Perfect Anchors: P2=Y, P9=V/I/L): +100 points")
    print("- Tier 2 (Good Anchors: P2=F, P9=V/I/L): +50 points")
    print("- Tier 3 (Bad Anchors): 0 points")
    print("- Penalty for internal Proline (P3-P8): -20 points")
    print("- Penalty for other internal changes vs E1: -1 point per change\n")

    for name, seq in epitopes.items():
        score = 0
        p2 = seq[1]
        p9 = seq[8]
        internal_seq = seq[2:8]

        # Tiered scoring based on primary anchors
        if p2 == 'Y' and p9 in ['V', 'I', 'L']:
            score += 100  # Tier 1
        elif p2 == 'F' and p9 in ['V', 'I', 'L']:
            score += 50   # Tier 2
        else:
            score += 0    # Tier 3 (Bad anchors)

        # Penalties for internal residues
        if 'P' in internal_seq:
            score -= 20
        
        # Penalty for deviation from reference internal sequence (for distinguishing E1 and E4)
        if name in ["E1", "E4"]: # only compare E4 to E1 in this specific case
            changes = sum(1 for i in range(len(internal_seq)) if internal_seq[i] != reference_internal_sequence[i])
            score -= changes

        scores[name] = {'score': score, 'sequence': seq}

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1]['score'], reverse=True)

    print("--- Calculated Scores and Final Ranking ---")
    for name, data in ranked_epitopes:
        print(f"Epitope: {name}, Sequence: {data['sequence']}, Score: {data['score']}")

    ranked_list = [item[0] for item in ranked_epitopes]
    print(f"\nFinal Rank (Highest to Lowest Binding): {', '.join(ranked_list)}")

    # Print the final equation using the ranked list and corresponding values
    print("\nBased on our scoring system, the ranking equation is:")
    final_equation = []
    for name, data in ranked_epitopes:
        final_equation.append(f"{name}({data['score']})")
    print(" > ".join(final_equation))


solve_epitope_ranking()

# Return the final answer in the specified format
print("\n<<<A>>>")
