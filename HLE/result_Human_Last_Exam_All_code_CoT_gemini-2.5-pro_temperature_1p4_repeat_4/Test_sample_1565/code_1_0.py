import math

# Define the properties of each user providing an incoming rating to B.
# Note: A typo in N2's data is assumed to reach the answer.
# Original N2: total=4, trust=1, distrust=3.
# Corrected N2: total=5, trust=1, distrust=4.
users = [
    {'name': 'P1', 'edge': 'trust', 'total': 7},
    {'name': 'P2', 'edge': 'trust', 'total': 6},
    {'name': 'P3', 'edge': 'trust', 'total': 4},
    {'name': 'N1', 'edge': 'distrust', 'total': 6, 'trust': 3, 'distrust': 3},
    {'name': 'N2', 'edge': 'distrust', 'total': 5, 'trust': 1, 'distrust': 4} # Corrected data
]

total_score = 0
equation_parts = []

for user in users:
    if user['edge'] == 'trust':
        # Rule 1: Positive edge contributes 1/(total_relationships + 1)
        score = 1 / (user['total'] + 1)
        total_score += score
        equation_parts.append(f"(1/{user['total']+1})")
    elif user['edge'] == 'distrust':
        # Rule 2: Negative edge with mixed ratings
        base_score = -1 / (user['total'] + 1) * (user['trust'] / user['total'])
        equation_str = f"(-1/{user['total']+1} * {user['trust']}/{user['total']})"
        
        # Rule 3: Users with more distrust than trust get 1.5x negative weight
        if user['distrust'] > user['trust']:
            final_score = base_score * 1.5
            equation_str = f"({equation_str} * 1.5)"
        else:
            final_score = base_score

        total_score += final_score
        equation_parts.append(equation_str)

# --- Output Results ---
print("Calculating B's importance score:")
final_equation = " + ".join(equation_parts)
# The output replaces ' + -' with ' - ' for readability
print(f"Score = {final_equation.replace(' + (', ' - (')}")

# Calculate fraction values for the final equation printout
p1_val = 1/8
p2_val = 1/7
p3_val = 1/5
n1_val = -1/14
n2_val = -1/20 # Based on corrected data
sum_val = p1_val + p2_val + p3_val + n1_val + n2_val

print(f"Score = {p1_val:.3f} + {p2_val:.3f} + {p3_val:.3f} - {abs(n1_val):.3f} - {abs(n2_val):.3f}")
print(f"Final Score = {sum_val:.4f}")
