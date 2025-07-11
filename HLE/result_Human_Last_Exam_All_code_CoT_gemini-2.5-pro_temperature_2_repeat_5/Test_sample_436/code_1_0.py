import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd
    by scoring them against a known binding motif.
    """
    # Define the epitopes
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }

    # Simplified scoring matrix based on H2-Kd binding motifs.
    # Positive values are favorable, negative values are unfavorable.
    # A value of 0 is assigned to neutral or unlisted residues.
    scoring_matrix = {
        2: {'Y': 3, 'F': 2},
        3: {'I': 2},
        4: {'P': -1},
        5: {'M': 2},
        9: {'V': 3, 'K': -5}
    }

    print("--- Epitope Binding Score Calculation for H2-Kd ---")
    
    # Calculate scores for each epitope
    results = {}
    for name, sequence in epitopes.items():
        score = 0
        calculation_str = []
        for i, amino_acid in enumerate(sequence):
            position = i + 1
            # Check if the position has specific scoring rules
            if position in scoring_matrix:
                # Get the score for the amino acid, default to 0 if not specified
                pos_score = scoring_matrix[position].get(amino_acid, 0)
                if pos_score != 0:
                    score += pos_score
                    calculation_str.append(f"P{position}({amino_acid}):{pos_score}")
        
        results[name] = score
        print(f"\n{name} ({sequence}):")
        print(f"  Calculation: {' + '.join(calculation_str).replace('+ -', '- ')}")
        print(f"  Final Score: {score}")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(results.items(), key=lambda item: item[1], reverse=True)

    print("\n--- Final Ranking (Highest to Lowest Binding Affinity) ---")
    ranked_list = [name for name, score in ranked_epitopes]
    print(" -> ".join(ranked_list))
    
rank_epitopes()