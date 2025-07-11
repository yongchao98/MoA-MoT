import math

# Step 1: Define the data for each user connecting to B
users = [
    # Trust Incoming (Positive Edges)
    # P1 (7 total): 4 trust-senior, 2 trust-peer -> 6 trust, 1 distrust
    {'name': 'P1', 'edge_type': 'positive', 'total_ratings': 7, 'trust_ratings': 6},
    # P2 (6 total): 3 trust-junior, 1 trust-senior -> 4 trust, 2 distrust
    {'name': 'P2', 'edge_type': 'positive', 'total_ratings': 6, 'trust_ratings': 4},
    # P3 (4 total): 2 trust-peer -> 2 trust, 2 distrust
    {'name': 'P3', 'edge_type': 'positive', 'total_ratings': 4, 'trust_ratings': 2},
    
    # Distrust Incoming (Negative Edges)
    # N1 (6 total): 3 trust-senior -> 3 trust, 3 distrust
    {'name': 'N1', 'edge_type': 'negative', 'total_ratings': 6, 'trust_ratings': 3},
    # N2 (4 total): 1 trust-peer -> 1 trust, 3 distrust
    {'name': 'N2', 'edge_type': 'negative', 'total_ratings': 4, 'trust_ratings': 1},
]

total_score = 0
contributions = []

# Steps 2 & 3: Calculate contribution for each user
for user in users:
    score = 0
    total_ratings = user['total_ratings']
    trust_ratings = user['trust_ratings']
    
    if user['edge_type'] == 'positive':
        # Rule 1: Positive edge contributes 1/(total_relationships + 1)
        score = 1 / (total_ratings + 1)
    else: # Negative edge
        # Rule 2: Negative edge with mixed ratings: -1 /(total relationships+1) x (trust_ratings/total)
        score = -1 / (total_ratings + 1) * (trust_ratings / total_ratings)
        
        # Rule 3: Users with more distrust than trust ratings get 1.5x negative weight
        distrust_ratings = total_ratings - trust_ratings
        if distrust_ratings > trust_ratings:
            score *= 1.5
            
    contributions.append(score)
    total_score += score

# Step 4 & 5: Print the final equation and score
print("B's FAT score is the sum of contributions from each connecting user:")
print("Final Equation:")

# Format the equation string with each number
equation_str = " + ".join([f"({c:.4f})" if c < 0 else f"{c:.4f}" for c in contributions])
print(f"{equation_str} = {total_score:.4f}")

print(f"\nFinal Calculated Score: {total_score:.4f}")

# Compare with the options provided
options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
closest_option = min(options, key=lambda k: abs(options[k] - total_score))
print(f"The calculated score is closest to option {closest_option}: {options[closest_option]}")