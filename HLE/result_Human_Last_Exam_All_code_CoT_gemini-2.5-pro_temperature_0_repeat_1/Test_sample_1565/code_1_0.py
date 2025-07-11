import collections

# Define a structure to hold user data
User = collections.namedtuple('User', ['name', 'type', 'total_ratings', 'trust_ratings', 'distrust_ratings'])

# Data for users connecting to B
users_data = [
    # Trust Incoming
    User(name='P1', type='trust', total_ratings=7, trust_ratings=6, distrust_ratings=1),
    User(name='P2', type='trust', total_ratings=6, trust_ratings=4, distrust_ratings=2),
    User(name='P3', type='trust', total_ratings=4, trust_ratings=2, distrust_ratings=2),
    # Distrust Incoming
    User(name='N1', type='distrust', total_ratings=6, trust_ratings=3, distrust_ratings=3),
    User(name='N2', type='distrust', total_ratings=4, trust_ratings=1, distrust_ratings=3),
]

total_score = 0
score_components = []

print("Calculating B's FAT score step-by-step:\n")

for user in users_data:
    score = 0
    if user.type == 'trust':
        # Rule 1: Positive edge contributes 1/(total_relationships + 1)
        score = 1 / (user.total_ratings + 1)
        print(f"Contribution from {user.name} (Trust): 1 / ({user.total_ratings} + 1) = {score}")
    elif user.type == 'distrust':
        # Rule 2: Negative edge with mixed ratings: -1 /(total relationships+1) x (trust_ratings/total)
        score = -1 / (user.total_ratings + 1) * (user.trust_ratings / user.total_ratings)
        
        # Rule 3: Users with more distrust than trust ratings get 1.5x negative weight
        if user.distrust_ratings > user.trust_ratings:
            print(f"Contribution from {user.name} (Distrust, Rule 3 applies): (-1 / ({user.total_ratings} + 1) * ({user.trust_ratings} / {user.total_ratings})) * 1.5 = {score * 1.5}")
            score *= 1.5
        else:
            print(f"Contribution from {user.name} (Distrust): -1 / ({user.total_ratings} + 1) * ({user.trust_ratings} / {user.total_ratings}) = {score}")

    score_components.append(score)
    total_score += score

print("\nB's total FAT score is the sum of all contributions:")
# Construct the final equation string, handling negative numbers
equation_parts = [f"({s})" if s < 0 else str(s) for s in score_components]
equation_str = " + ".join(equation_parts)
print(f"Final Score = {equation_str}")

print(f"Final Score = {total_score}")

# Compare with options
options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
closest_option = min(options, key=lambda k: abs(options[k] - total_score))

print(f"\nThe calculated score {total_score:.3f} is closest to option {closest_option}) {options[closest_option]}.")
print(f"\n{closest_option}")
<<<A>>>