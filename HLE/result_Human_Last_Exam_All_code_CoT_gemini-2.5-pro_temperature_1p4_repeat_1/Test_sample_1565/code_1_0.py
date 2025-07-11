import math

# Data for incoming connections to user B
# User: (is_trust_edge, total_ratings, trust_ratings, distrust_ratings)
connections = {
    'P1': (True, 7, 4 + 2, 1),
    'P2': (True, 6, 3 + 1, 2),
    'P3': (True, 4, 2, 1 + 1),
    'N1': (False, 6, 3, 2 + 1),
    'N2': (False, 4, 1, 3),
}

total_score = 0
equation_parts = []
print("B's importance score calculation (FAT):")
print("-" * 35)

# Calculate contribution for each connection
for user, data in connections.items():
    is_trust_edge, total_ratings, trust_ratings, distrust_ratings = data
    contribution = 0
    
    if is_trust_edge:
        # Rule 1: Positive edge
        contribution = 1 / (total_ratings + 1)
        print(f"{user} (Trust):  1 / ({total_ratings} + 1) = {contribution:.4f}")
        equation_parts.append(f"1/{total_ratings+1}")
    else:
        # Rule 2: Negative edge with mixed ratings
        # All negative edges in this problem have mixed ratings
        contribution = -1 / (total_ratings + 1) * (trust_ratings / total_ratings)
        base_calculation_str = f"-1 / ({total_ratings} + 1) * ({trust_ratings} / {total_ratings})"
        
        # Rule 3: Check for 1.5x negative weight
        multiplier = 1.0
        if distrust_ratings > trust_ratings:
            multiplier = 1.5
            contribution *= multiplier
            print(f"{user} (Distrust): {multiplier} * [{base_calculation_str}] = {contribution:.4f}")
            equation_parts.append(f"{multiplier}*(-1/({(total_ratings+1)*total_ratings/trust_ratings}))")

        else:
            print(f"{user} (Distrust): {base_calculation_str} = {contribution:.4f}")
            equation_parts.append(f"-1/({(total_ratings+1)*total_ratings/trust_ratings})")

    total_score += contribution

# Print the final summary
print("-" * 35)
print("Final Equation:")
# The equation string composition is simplified for readability
print(f"Score = (1/8) + (1/7) + (1/5) - (1/14) - (1.5/20)")

print("\nFinal Calculation:")
calculation_str = " + ".join(f"{p:.4f}" for p in [
    1/8, 1/7, 1/5, -1/14, -1.5/20
])
print(f"Score = {calculation_str} = {total_score:.4f}")

print("\nThe calculated score is approximately 0.3214.")
print("The closest option is A) 0.35.")