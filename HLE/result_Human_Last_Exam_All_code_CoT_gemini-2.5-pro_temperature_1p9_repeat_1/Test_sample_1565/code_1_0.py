import math

def calculate_fat_score():
    """
    Calculates and prints B's FAT importance score based on the provided rules.
    """
    # Data for users with relationships to B
    # Format: {'name', 'is_positive_edge', 'total_rels', 'trust_rels', 'distrust_rels'}
    users = [
        {'name': 'P1', 'is_positive': True,  'total': 7, 'trust': 6, 'distrust': 1},
        {'name': 'P2', 'is_positive': True,  'total': 6, 'trust': 4, 'distrust': 2},
        {'name': 'P3', 'is_positive': True,  'total': 4, 'trust': 2, 'distrust': 2},
        {'name': 'N1', 'is_positive': False, 'total': 6, 'trust': 3, 'distrust': 3},
        {'name': 'N2', 'is_positive': False, 'total': 4, 'trust': 1, 'distrust': 3}
    ]

    total_score = 0
    scores = []

    for user in users:
        total_rels = user['total']
        if user['is_positive']:
            # Rule 1: Positive edge contributes 1/(total_relationships + 1)
            score = 1 / (total_rels + 1)
        else:
            # Negative edge
            trust_ratings = user['trust']
            distrust_ratings = user['distrust']
            
            # Rule 2: Negative edge with mixed ratings
            # -1 /(total relationships+1) x (trust_ratings/total)
            score = -1 / (total_rels + 1) * (trust_ratings / total_rels)
            
            # Rule 3: Users with more distrust than trust get 1.5x negative weight
            if distrust_ratings > trust_ratings:
                score *= 1.5
        
        scores.append(score)
        total_score += score
    
    # Output the required final equation
    print("B's importance score equation:")
    equation_parts = [f"({s:.3f})" if s < 0 else f"{s:.3f}" for s in scores]
    equation_str = " + ".join(equation_parts)
    print(f"Score = {equation_str}")
    
    # Print the final result
    print(f"\nFinal Calculated Score: {total_score}")

    # Determine closest option
    options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
    closest_option = min(options.items(), key=lambda item: abs(item[1] - total_score))
    print(f"The calculated score {total_score:.3f} is closest to option {closest_option[0]}: {closest_option[1]}")

calculate_fat_score()