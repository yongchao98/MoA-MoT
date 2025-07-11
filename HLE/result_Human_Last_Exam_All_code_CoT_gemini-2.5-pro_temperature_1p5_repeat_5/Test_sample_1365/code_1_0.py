# Parameters for the problem, you can change these values
n = 10  # Total number of experts
c = 5   # Mistake threshold

# --- Calculation of the Upper Bound ---
# Based on the derived formula: Upper Bound = n * c - 1

print(f"Calculating the upper bound for n={n} and c={c}:")
print("--------------------------------------------------")

# Part 1: Bound from the true expert's mistakes.
# The true expert makes strictly fewer than c mistakes.
max_true_expert_mistakes = c - 1
print(f"The true expert makes at most c - 1 = {max_true_expert_mistakes} mistakes.")

# Part 2: Bound from the false experts' mistakes.
# There are n-1 false experts, each making at most c mistakes.
num_false_experts = n - 1
max_total_false_mistakes = num_false_experts * c
print(f"The {num_false_experts} false experts make at most (n - 1) * c = {max_total_false_mistakes} mistakes in total.")

# Part 3: Combine the bounds for the final result.
# The algorithm's mistakes M are bounded by the sum of these two values.
upper_bound = max_total_false_mistakes + max_true_expert_mistakes
print("\nThe algorithm's total mistakes (M) is bounded by the sum:")
print(f"M <= (Total false expert mistakes) + (True expert mistakes)")
print(f"M <= {max_total_false_mistakes} + {max_true_expert_mistakes}")
print(f"M <= {upper_bound}")

print("\n--- Final Equation ---")
print("The simplified final equation for the upper bound is: M <= n * c - 1")
print(f"Plugging in the numbers: M <= {n} * {c} - 1 = {upper_bound}")
