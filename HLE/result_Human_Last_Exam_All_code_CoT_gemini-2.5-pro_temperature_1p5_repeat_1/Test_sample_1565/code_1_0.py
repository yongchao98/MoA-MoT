# User profiles: [total_ratings, trust_ratings, distrust_ratings]
users = {
    'P1': {'type': 'trust', 'profile': [7, 6, 1]},
    'P2': {'type': 'trust', 'profile': [6, 4, 2]},
    'P3': {'type': 'trust', 'profile': [4, 2, 2]},
    'N1': {'type': 'distrust', 'profile': [6, 3, 3]},
    'N2': {'type': 'distrust', 'profile': [4, 1, 3]}
}

total_score = 0
equation_terms = []

# Calculate score for each user
for user, data in users.items():
    total_rels = data['profile'][0]
    trust_ratings = data['profile'][1]
    distrust_ratings = data['profile'][2]
    
    score = 0
    term_str = ""

    if data['type'] == 'trust':
        # Rule 1: Positive edge
        score = 1 / (total_rels + 1)
        term_str = f"(1/{total_rels + 1})"
    else: # type == 'distrust'
        # Rule 2: Negative edge with mixed ratings
        # All negative sources (N1, N2) have mixed ratings
        base_score = -1 / (total_rels + 1) * (trust_ratings / total_rels)
        score = base_score
        
        # Rule 3: Extra weight for more distrust than trust
        if distrust_ratings > trust_ratings:
            score *= 1.5
            # For pretty printing the fraction
            if user == 'N2': # Specifically for N2 in this problem
                term_str = f"(-3/40)"
        else:
            # For pretty printing N1's fraction
             if user == 'N1':
                term_str = f"(-1/14)"

    total_score += score
    equation_terms.append(term_str)

# Final output
equation = " + ".join(equation_terms).replace("+ -", "- ")
print(f"B's importance score calculation:")
# Print the equation with fractions for clarity, as derived in the thinking steps
print("Score = 1/8 + 1/7 + 1/5 - 1/14 - 3/40")
# Print the final result
print(f"Score = {1/8:.4f} + {1/7:.4f} + {1/5:.4f} - {1/14:.4f} - {3/40:.4f}")
print(f"Total Score = {total_score:.4f}")

# Based on the calculation (0.3214), the closest option is A (0.35)
print("\nThe calculated score is approximately 0.3214. The closest option is A.")