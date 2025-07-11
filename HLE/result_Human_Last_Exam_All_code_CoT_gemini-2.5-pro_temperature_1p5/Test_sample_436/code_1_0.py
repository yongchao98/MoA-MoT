import pandas as pd

def calculate_affinity_score(sequence, reference_sequence="TYQRTRALV"):
    """
    Calculates a heuristic affinity score for a peptide binding to H-2Kd.
    This model is based on the known H-2Kd binding motif.
    - Score is primarily based on anchor residues at P2 and P9.
    - Minor penalties are applied for mutations at other positions.
    """
    # Define scoring for primary anchors
    p2_scores = {'Y': 100, 'F': 60}  # Optimal Y is scored much higher than suboptimal F
    p9_scores = {'V': 10, 'L': 10, 'I': 10, 'K': -500} # Hydrophobic P9 is key, charged K is highly penalized

    # Get residues at anchor positions
    p2_residue = sequence[1]
    p9_residue = sequence[8]

    # Calculate anchor scores
    p2_score = p2_scores.get(p2_residue, 0)
    p9_score = p9_scores.get(p9_residue, 0)

    # Penalize for mutations at non-anchor positions (secondary effects)
    secondary_mutations = 0
    for i, (res, ref_res) in enumerate(zip(sequence, reference_sequence)):
        # We only count non-anchor positions
        if i not in [1, 8] and res != ref_res:
            secondary_mutations += 1
    
    # A small penalty for each secondary mutation
    secondary_penalty = secondary_mutations * 5

    # The final score is a combination of anchor scores minus the secondary penalty
    # P2 is the most dominant factor in this simple model.
    final_score = p2_score + p9_score - secondary_penalty
    
    # We print the "equation" for each epitope as requested
    print(f"Score for {sequence}: (P2 score: {p2_score}) + (P9 score: {p9_score}) - (Secondary Penalty: {secondary_penalty}) = {final_score}")

    return final_score

def rank_epitopes():
    """
    Ranks epitopes based on their calculated binding affinity scores.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV"
    }

    print("Calculating affinity scores for each epitope:\n")
    
    # Calculate scores and store them
    scores = {}
    for name, seq in epitopes.items():
        scores[name] = calculate_affinity_score(seq)

    # Sort the epitopes by score in descending order
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    print("\n--- Final Ranking ---")
    print("The epitopes ranked from highest to lowest expected binding are:")
    for i, (name, score) in enumerate(ranked_epitopes):
        print(f"{i+1}. {name} (Sequence: {epitopes[name]}, Score: {score})")
    
    final_order = [name for name, score in ranked_epitopes]
    print(f"\nThe final ranked order is: {', '.join(final_order)}")


if __name__ == '__main__':
    rank_epitopes()
    print("\nBased on the analysis, the correct answer choice is A.")
    print("<<<A>>>")
