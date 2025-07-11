# Step 1: Define properties of the connecting users
# B's incoming trust relationships
positive_sources = ["P1", "P2", "P3"]

# B's incoming distrust relationships
# N1: 3 trust-senior, 2 distrust-peer, 1 distrust-junior -> 3 trust, 3 distrust, 6 total
n1_props = {'trust': 3, 'distrust': 3, 'total': 6}
# N2: 1 trust-peer, 3 distrust-junior -> 1 trust, 3 distrust, 4 total
n2_props = {'trust': 1, 'distrust': 3, 'total': 4}

# Step 2: Calculate B's total in-degree
b_in_degree = len(positive_sources) + 2 # P1,P2,P3 and N1,N2
denominator = b_in_degree + 1

# Step 3: Calculate contribution from positive edges (Rule 1)
p_contribution = 1 / denominator

# Step 4: Calculate contribution from negative edges (Rules 2 & 3)

# For N1
n1_base_contribution = -1 / denominator * (n1_props['trust'] / n1_props['total'])
# Rule 3 check for N1: distrust (3) is not greater than trust (3)
n1_final_contribution = n1_base_contribution

# For N2
n2_base_contribution = -1 / denominator * (n2_props['trust'] / n2_props['total'])
# Rule 3 check for N2: distrust (3) is greater than trust (1), apply 1.5x weight
n2_final_contribution = n2_base_contribution * 1.5

# Step 5: Sum all contributions
total_score = (len(positive_sources) * p_contribution) + n1_final_contribution + n2_final_contribution

# Step 6: Print the final equation and result
print("B's Importance Score Calculation:")
print(f"B's total in-degree = {len(positive_sources)} (trust) + 2 (distrust) = {b_in_degree}")
print(f"Denominator used for all calculations = (in-degree + 1) = {b_in_degree + 1}")
print("-" * 30)

# Print equation components
p1_str = f"P1: (1/{denominator})"
p2_str = f"P2: (1/{denominator})"
p3_str = f"P3: (1/{denominator})"
n1_str = f"N1: (-1/{denominator} * {n1_props['trust']}/{n1_props['total']})"
n2_str = f"N2: (-1/{denominator} * {n2_props['trust']}/{n2_props['total']}) * 1.5"

print(f"Final Score = {p1_str} + {p2_str} + {p3_str} + {n1_str} + {n2_str}")
print(f"Final Score = {p_contribution:.4f} + {p_contribution:.4f} + {p_contribution:.4f} + {n1_final_contribution:.4f} + {n2_final_contribution:.4f}")
print(f"Final Score = {total_score:.4f}")
print("-" * 30)
print("The calculated score is approximately 0.354, which corresponds to option A.")
