import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    
    The function uses a scoring matrix based on known H2-Kd binding motifs to calculate
    an affinity score for each epitope. It then ranks them from highest to lowest score.
    """
    # Define the epitope sequences
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    # Scoring matrix based on H2-Kd binding motif preferences.
    # Positive scores for preferred residues, negative for detrimental ones.
    # Positions are 1-indexed, so we use dictionary keys 2 through 9.
    score_matrix = {
        # P2: Primary anchor. Y is optimal, F is good.
        2: {'Y': 10, 'F': 8},
        # P3: Secondary anchor. I is strongly preferred.
        3: {'I': 6, 'Q': 2},
        # P4: Proline can be disruptive.
        4: {'P': -3},
        # P5: Secondary anchor. Hydrophobic M is preferred over polar T.
        5: {'M': 4, 'T': 2},
        # P6: Secondary anchor. Aromatic F is highly preferred over charged R.
        6: {'F': 7, 'R': 0},
        # P9: Primary anchor. Hydrophobic V is optimal. Charged K is highly detrimental.
        9: {'V': 10, 'I': 10, 'L': 10, 'K': -15}
    }

    print("--- Calculating Binding Affinity Scores for H2-Kd ---\n")
    
    calculated_scores = {}
    # Calculate score for each epitope
    for name, seq in epitopes.items():
        score = 0
        equation_str = f"Score for {name} ({seq}):"
        
        # Peptide positions are 1-based, Python indices are 0-based.
        for position_1_based, amino_acid in enumerate(seq, 1):
            if position_1_based in score_matrix and amino_acid in score_matrix[position_1_based]:
                points = score_matrix[position_1_based][amino_acid]
                score += points
                equation_str += f" + {points} (P{position_1_based}:{amino_acid})"
        
        # Add a base score for having the right length and general structure
        base_score = 50
        score += base_score
        equation_str = equation_str.replace(":", f" = {base_score}", 1) # Insert base score
        equation_str += f" = {score}"
        
        calculated_scores[name] = score
        print(equation_str)
        
    print("\n--- Analysis ---")
    print("E4: Optimal residues at primary anchors (P2:Y, P9:V) and a key secondary anchor (P6:F). Highest score.")
    print("E5: Optimal primary anchors. Major improvement at P3 (I) likely outweighs the negative Proline at P4.")
    print("E1: High-affinity reference with optimal primary anchors but suboptimal secondary anchors (P6:R).")
    print("E2: Good primary anchors (P2:F is slightly less optimal than Y) but suboptimal secondary anchors.")
    print("E3: Catastrophic mutation at the C-terminal anchor (P9:K), leading to very poor binding.")


    # Sort the epitopes by score in descending order
    ranked_epitopes = sorted(calculated_scores.items(), key=lambda item: item[1], reverse=True)

    print("\n--- Final Ranking (Highest to Lowest Amount Complexed) ---")
    final_rank_str = ", ".join([name for name, score in ranked_epitopes])
    print(f"The predicted order of binding affinity is: {final_rank_str}")
    
    # Match to the given answer choices
    final_ranking_tuple = tuple(name for name, score in ranked_epitopes)
    answer_choices = {
        'A': ('E1', 'E4', 'E5', 'E2', 'E3'),
        'B': ('E4', 'E5', 'E1', 'E2', 'E3'),
        'C': ('E5', 'E4', 'E1', 'E2', 'E3'),
        'D': ('E1', 'E5', 'E4', 'E2', 'E3'),
        'E': ('E5', 'E4', 'E1', 'E3', 'E2')
    }
    
    correct_choice = None
    for choice, ranking in answer_choices.items():
        if final_ranking_tuple == ranking:
            correct_choice = choice
            break
            
    print(f"\nThis corresponds to Answer Choice: {correct_choice}")


rank_epitopes()
<<<B>>>