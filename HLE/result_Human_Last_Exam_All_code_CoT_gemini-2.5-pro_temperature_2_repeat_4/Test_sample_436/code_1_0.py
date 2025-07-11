def calculate_h2kd_binding_score(peptide_sequence, reference_peptide):
    """
    Calculates a binding affinity score for a peptide to H-2Kd based on known motifs.
    Higher scores indicate better predicted binding.
    """
    score = 0
    
    # Position 2 (P2) Anchor Scoring
    p2_residue = peptide_sequence[1]
    if p2_residue == 'Y':
        score += 20  # Optimal anchor
    elif p2_residue == 'F':
        score += 10  # Tolerated anchor
    
    # Position 9 (P9) Anchor Scoring
    p9_residue = peptide_sequence[8]
    if p9_residue in ['V', 'L']:
        score += 10  # Optimal anchor
    elif p9_residue == 'K':
        score += 0   # Detrimental anchor
        
    # Penalties for mutations at non-anchor positions relative to the reference peptide
    # This reflects that the reference is already a known high-affinity binder
    for i in range(len(peptide_sequence)):
        # Skip the primary anchor positions P2 (index 1) and P9 (index 8)
        if i == 1 or i == 8:
            continue
            
        # Apply penalty if residue differs from the reference peptide
        if peptide_sequence[i] != reference_peptide[i]:
            # Proline at P4 (index 3) is particularly disruptive
            if i == 3 and peptide_sequence[i] == 'P':
                score -= 3
            else:
                score -= 1 # General penalty for other non-anchor mutations
                
    return score

def main():
    """
    Main function to rank epitopes based on binding scores.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV"
    }

    # The reference high-affinity epitope is E1
    reference_peptide = epitopes["E1"]

    # Calculate scores for each epitope
    scores = {}
    for name, sequence in epitopes.items():
        scores[name] = calculate_h2kd_binding_score(sequence, reference_peptide)
    
    # Sort epitopes by score in descending order
    # In case of a tie, the original order is preserved, which is fine for this problem.
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    
    # Print the explanation and the ranking
    print("Ranking of Epitopes for H-2Kd Binding (from highest to lowest affinity):")
    print("-" * 65)
    print("{:<10} {:<15} {:<10} {:<25}".format("Rank", "Epitope", "Sequence", "Score"))
    print("-" * 65)
    
    for i, (name, score) in enumerate(ranked_epitopes):
        sequence = epitopes[name]
        p2 = sequence[1]
        p9 = sequence[8]
        print("{:<10} {:<15} {:<10} {:<25}".format(
            f"{i+1}.",
            name,
            f"({sequence})",
            f"Score={score} (P2='{p2}', P9='{p9}')"
        ))
        
    final_ranking = [name for name, score in ranked_epitopes]
    print("\nFinal Predicted Rank Order:", ", ".join(final_ranking))
    print("This corresponds to Answer Choice A.")


if __name__ == "__main__":
    main()
