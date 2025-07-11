import pandas as pd

def count_internal_differences(seq1, seq2):
    """Counts differences in the internal parts of two 9-mer sequences (P3-P8)."""
    internal1 = seq1[2:8]
    internal2 = seq2[2:8]
    return sum(1 for a, b in zip(internal1, internal2) if a != b)

def rank_epitopes():
    """
    Ranks MHC binding epitopes based on a scoring model derived from immunological principles.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    reference_epitope = 'TYQRTRALV'
    results = []

    print("Ranking epitopes by predicted binding affinity to H2-Kd...\n")
    print("Scoring Model:")
    print("- Base Score: +100")
    print("- P2='Y' (Optimal): +100")
    print("- P2='F' (Sub-optimal): +0")
    print("- P9='V' (Optimal): +50")
    print("- P9='K' (Deleterious): -500")
    print("- Internal Proline: -20")
    print("- Each Internal Amino Acid Difference from Reference E1: -10\n")


    for name, seq in epitopes.items():
        score = 100
        p2_score = 0
        p9_score = 0
        proline_penalty = 0
        diff_penalty = 0

        # P2 Anchor Score
        if seq[1] == 'Y':
            p2_score = 100
        elif seq[1] == 'F':
            p2_score = 0

        # P9 Anchor Score
        if seq[8] == 'V':
            p9_score = 50
        elif seq[8] == 'K':
            p9_score = -500
        
        # Internal Proline Penalty
        if 'P' in seq[2:8]:
            proline_penalty = -20

        # Internal Differences Penalty
        diff_count = count_internal_differences(seq, reference_epitope)
        diff_penalty = -10 * diff_count
        
        total_score = score + p2_score + p9_score + proline_penalty + diff_penalty
        
        results.append({
            'Epitope': name,
            'Sequence': seq,
            'Score': total_score,
            'Calculation': f"{score} (Base) + {p2_score} (P2) + {p9_score} (P9) + {proline_penalty} (Proline) + {diff_penalty} (Diffs)"
        })

    # Sort epitopes by score in descending order
    ranked_results = sorted(results, key=lambda x: x['Score'], reverse=True)

    # Print the ranked list with calculations
    for i, result in enumerate(ranked_results):
        print(f"{i+1}. {result['Epitope']} ({result['Sequence']}): Score = {result['Score']}")
        print(f"   Calculation: {result['Calculation']} = {result['Score']}\n")

    # Print the final ranked order as a simple string
    final_order = ", ".join([r['Epitope'] for r in ranked_results])
    print(f"Final Ranked Order: {final_order}")

if __name__ == '__main__':
    rank_epitopes()
