import collections

# Define the properties of each user connecting to B.
# 'total_out_ratings' is the user's total number of ratings given.
# 'trust_out_ratings' is the number of trust ratings given by the user.
User = collections.namedtuple('User', ['name', 'edge_type', 'total_out_ratings', 'trust_out_ratings', 'distrust_out_ratings'])

users = [
    User(name='P1', edge_type='trust', total_out_ratings=7, trust_out_ratings=6, distrust_out_ratings=1),
    User(name='P2', edge_type='trust', total_out_ratings=6, trust_out_ratings=4, distrust_out_ratings=2),
    User(name='P3', edge_type='trust', total_out_ratings=4, trust_out_ratings=2, distrust_out_ratings=2),
    User(name='N1', edge_type='distrust', total_out_ratings=6, trust_out_ratings=3, distrust_out_ratings=3),
    User(name='N2', edge_type='distrust', total_out_ratings=4, trust_out_ratings=1, distrust_out_ratings=3),
]

# The number of incoming relationships to B (B's in-degree)
b_in_degree = len(users)

total_score = 0
equation_parts = []

print("Calculating B's FAT Score:")
print(f"B's total incoming relationships = {b_in_degree}\n")

for user in users:
    contribution = 0
    # Rule 1: Positive edge
    if user.edge_type == 'trust':
        contribution = 1 / (b_in_degree + 1)
        equation_parts.append(f"{contribution:.4f}")
        print(f"Contribution from {user.name} (Positive Edge):")
        print(f"  Formula: 1 / (total_relationships + 1)")
        print(f"  Calculation: 1 / ({b_in_degree} + 1) = {contribution:.4f}\n")

    # Rule 2 & 3: Negative edge
    elif user.edge_type == 'distrust':
        base_negative_score = -1 / (b_in_degree + 1) * (user.trust_out_ratings / user.total_out_ratings)
        
        # Rule 3: Check for 1.5x multiplier
        if user.distrust_out_ratings > user.trust_out_ratings:
            contribution = 1.5 * base_negative_score
            print(f"Contribution from {user.name} (Negative Edge with > Distrust):")
            print(f"  Formula: 1.5 * [-1 / (total_relationships + 1) * (trust_ratings / total)]")
            print(f"  Calculation: 1.5 * [-1 / ({b_in_degree} + 1) * ({user.trust_out_ratings} / {user.total_out_ratings})] = {contribution:.4f}\n")
        else:
            contribution = base_negative_score
            print(f"Contribution from {user.name} (Negative Edge):")
            print(f"  Formula: -1 / (total_relationships + 1) * (trust_ratings / total)")
            print(f"  Calculation: -1 / ({b_in_degree} + 1) * ({user.trust_out_ratings} / {user.total_out_ratings}) = {contribution:.4f}\n")

        equation_parts.append(f"({contribution:.4f})")
    
    total_score += contribution

# Print the final equation and result
equation_str = " + ".join(equation_parts).replace("+ (", "- ")
print("---")
print("Final Score Calculation:")
print(f"Equation: {equation_str} = {total_score:.4f}")
print(f"The calculated score is approximately {total_score:.2f}, which corresponds to option A.")