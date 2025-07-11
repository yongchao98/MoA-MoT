# User data provided in the problem
users = {
    'P1': {'type': 'trust', 'total_ratings': 7, 'trust_ratings': 6, 'distrust_ratings': 1},
    'P2': {'type': 'trust', 'total_ratings': 6, 'trust_ratings': 4, 'distrust_ratings': 2},
    'P3': {'type': 'trust', 'total_ratings': 4, 'trust_ratings': 2, 'distrust_ratings': 2},
    'N1': {'type': 'distrust', 'total_ratings': 6, 'trust_ratings': 3, 'distrust_ratings': 3},
    'N2': {'type': 'distrust', 'total_ratings': 4, 'trust_ratings': 1, 'distrust_ratings': 3}
}

total_score = 0
contributions = {}

# Calculate contributions for each user
for name, data in users.items():
    total_rels = data['total_ratings']
    
    # Rule 1: Positive edge contribution
    if data['type'] == 'trust':
        score = 1 / (total_rels + 1)
        contributions[name] = score
    
    # Rule 2 & 3: Negative edge contribution
    elif data['type'] == 'distrust':
        trust_ratings = data['trust_ratings']
        distrust_ratings = data['distrust_ratings']
        
        # This interpretation assumes 'total' in the fraction also uses '+1' for consistency
        # as it leads to one of the provided answers.
        # Formula: -1/(total+1) * (trust_ratings/(total+1))
        score = -1 / (total_rels + 1) * (trust_ratings / (total_rels + 1))
        
        # Rule 3: 1.5x negative weight for users with more distrust than trust
        if distrust_ratings > trust_ratings:
            score *= 1.5
        
        contributions[name] = score

total_score = sum(contributions.values())

# Build and print the final equation
equation_parts = []
for name, score in contributions.items():
    equation_parts.append(f"{score:.4f}")

equation_str = " + ".join(equation_parts).replace('+ -', '- ')
print("B's Importance Score Calculation:")
print(f"P1 Contribution (Trust): 1 / (7 + 1) = {contributions['P1']:.4f}")
print(f"P2 Contribution (Trust): 1 / (6 + 1) = {contributions['P2']:.4f}")
print(f"P3 Contribution (Trust): 1 / (4 + 1) = {contributions['P3']:.4f}")
print(f"N1 Contribution (Distrust): -1 / (6 + 1) * (3 / (6 + 1)) = {contributions['N1']:.4f}")
print(f"N2 Contribution (Distrust): [-1 / (4 + 1) * (1 / (4 + 1))] * 1.5 = {contributions['N2']:.4f}")
print("\nFinal Equation:")
print(f"{equation_str} = {total_score:.4f}")

# Find the closest option
options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
closest_option = min(options, key=lambda k: abs(options[k] - total_score))
print(f"\nThe calculated score {total_score:.4f} is closest to option {closest_option}: {options[closest_option]}")
<<<A>>>