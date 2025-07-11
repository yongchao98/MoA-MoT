import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.

    The scoring is based on established H2-Kd binding motifs:
    - P2 Anchor (Position 2): Tyrosine (Y) is optimal. Phenylalanine (F) is suboptimal.
    - P9 Anchor (Position 9): Valine (V) or Leucine (L) are optimal. Lysine (K) is highly disfavored.
    - Core Residues: Proline (P) in the core (P3-P8) can be disruptive.
    - The reference epitope E1 is used as a tie-breaker.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }
    reference_epitope_seq = 'TYQRTRALV'
    
    scores = {}
    
    print("Calculating binding scores based on H2-Kd motif rules...\n")
    
    for name, seq in epitopes.items():
        score = 0
        p2 = seq[1]
        p9 = seq[8]
        core = seq[2:8]
        
        # P2 anchor score
        p2_score = 0
        if p2 == 'Y':
            p2_score = 100
        elif p2 == 'F':
            p2_score = 50
        
        # P9 anchor score
        p9_score = 0
        if p9 in ['V', 'L']:
            p9_score = 100
        elif p9 == 'K':
            p9_score = -100 # Strong penalty for disfavored residue
            
        # Core disruption penalty
        proline_penalty = 0
        if 'P' in core:
            proline_penalty = -40 # Penalty for disruptive proline
        
        # Reference bonus for tie-breaking
        reference_bonus = 0
        if seq == reference_epitope_seq:
            reference_bonus = 1
            
        total_score = p2_score + p9_score + proline_penalty + reference_bonus
        scores[name] = total_score
        
        print(f"Epitope: {name} ({seq})")
        print(f"  P2 ('{p2}'): {p2_score} | P9 ('{p9}'): {p9_score} | Proline Penalty: {proline_penalty} | Ref Bonus: {reference_bonus}")
        print(f"  Total Score = {p2_score} + {p9_score} + {proline_penalty} + {reference_bonus} = {total_score}\n")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    
    print("Final Ranking (from highest to lowest amount complexed with H2-Kd):")
    final_rank_list = [name for name, score in ranked_epitopes]
    print(" -> ".join(final_rank_list))

rank_epitopes()
<<<A>>>