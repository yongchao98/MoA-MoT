# User data based on the problem description
# P = Positive edge (Trust Incoming), N = Negative edge (Distrust Incoming)
# (total_ratings, trust_ratings, distrust_ratings)
users = {
    "P1": (7, 6, 1),
    "P2": (6, 4, 2),
    "P3": (4, 2, 2),
    "N1": (6, 3, 3),
    "N2": (4, 1, 3)
}

# --- Calculation ---
contributions = {}
equation_parts = []

# Positive contributions (Rule 1)
for user_id in ["P1", "P2", "P3"]:
    total_ratings, _, _ = users[user_id]
    score = 1 / (total_ratings + 1)
    contributions[user_id] = score
    equation_parts.append(f"1/({total_ratings}+1)")

# Negative contributions (Rules 2 & 3)
for user_id in ["N1", "N2"]:
    total_ratings, trust_ratings, distrust_ratings = users[user_id]
    
    # Rule 2: Base weight for negative edge with mixed ratings
    base_weight = -1 / (total_ratings + 1) * (trust_ratings / total_ratings)
    
    equation_part = f"(-1/({total_ratings}+1) * {trust_ratings}/{total_ratings}"

    # Rule 3: 1.5x multiplier for more distrust than trust
    if distrust_ratings > trust_ratings:
        final_weight = base_weight * 1.5
        equation_part += " * 1.5)"
    else:
        final_weight = base_weight
        equation_part += ")"

    contributions[user_id] = final_weight
    equation_parts.append(equation_part)

# Summing up all contributions
total_score = sum(contributions.values())

# --- Output ---
print("Calculating B's importance score (FAT):")
print("\nIndividual Contributions:")
for user, score in contributions.items():
    print(f"{user}: {score:.4f}")

print("\nFinal Equation:")
equation_string = " + ".join(equation_parts)
print(f"Score = {equation_string}")

# Print equation with numeric values
numeric_parts = [f"{score:.4f}" for score in contributions.values()]
print(f"Score = {' + '.join(numeric_parts).replace('+ -', '- ')}")


print(f"\nFinal Calculated Score: {total_score:.4f}")

# The calculated score 0.3214 is closest to option A (0.35)
print("\nThe closest answer option is A.")