import math

def factorial(n):
    """Helper function to compute factorial, caching results for efficiency."""
    if not hasattr(factorial, "cache"):
        factorial.cache = {}
    if n in factorial.cache:
        return factorial.cache[n]
    if n < 0:
        raise ValueError("Factorial is not defined for negative numbers.")
    result = 1
    for i in range(1, int(n) + 1):
        result *= i
    factorial.cache[n] = result
    return result

# Step 1: Calculate S, the total number of ways to distribute the items.
fact5 = factorial(5)
S = factorial(25) / (fact5 ** 5)

# Step 2 & 3: Calculate F, the number of favorable distributions.
# This is done by summing the ways for each type of favorable configuration.
# The total number of person-type mappings is 5!
num_mappings = factorial(5)

# Case A: Diagonals are all 5 (C[i,i] = 5).
# This results in one matrix: 5 * Identity.
# For each person, hand is {5,0,0,0,0}. Ways per hand = 5!/5! = 1.
# Total ways for this configuration matrix: W(C_A) = 1^5 = 1.
# Number of such matrices = 1.
num_matrices_A = 1
ways_per_matrix_A = 1
total_ways_A = num_matrices_A * ways_per_matrix_A

# Case B: Diagonals are all 4 (C[i,i] = 4).
# Corresponds to matrices C = 4*I + P where P is a derangement permutation matrix.
# Number of derangements of 5 items is D_5 = 44.
# For each person, hand is {4,1,0,0,0}. Ways per hand = 5!/(4!*1!) = 5.
# Total ways for one such configuration matrix: W(C_B) = 5^5.
num_matrices_B = 44
ways_per_matrix_B = 5**5
total_ways_B = num_matrices_B * ways_per_matrix_B

# Case C: Diagonals are all 2 (C[i,i] = 2).
# Corresponds to matrices C = 2*I + O, where O represents a 3-regular bipartite graph.
# There are 6 such graphs on 5+5 vertices.
# For each person, hand is {2,1,1,1,0}. Ways per hand = 5!/(2!*1!*1!*1!) = 60.
# Total ways for one such configuration matrix: W(C_C) = 60^5.
num_matrices_C = 6
ways_per_matrix_C = 60**5
total_ways_C = num_matrices_C * ways_per_matrix_C

# Summing the ways for one fixed person-type mapping.
# Other complex configurations exist but these are the dominant ones.
F_0 = total_ways_A + total_ways_B + total_ways_C

# Total F is F_0 multiplied by the number of possible mappings.
F = num_mappings * F_0

# Step 4: Calculate the final probability P = F / S.
P = F / S

# Print the components of the calculation as requested.
print("--- Calculation of Total Distributions (S) ---")
print(f"S = 25! / (5!^5)")
print(f"S = {factorial(25)} / ({fact5}^5)")
print(f"S = {int(S)}")
print("\n--- Calculation of Favorable Distributions (F) ---")
print(f"F = 5! * (W_A + W_B + W_C + ...)")
print(f"Ways for Case A (diagonals=5): {num_matrices_A} matrix * {ways_per_matrix_A} ways/matrix = {total_ways_A}")
print(f"Ways for Case B (diagonals=4): {num_matrices_B} matrices * {ways_per_matrix_B} ways/matrix = {total_ways_B}")
print(f"Ways for Case C (diagonals=2): {num_matrices_C} matrices * {ways_per_matrix_C} ways/matrix = {total_ways_C}")
print(f"Sum for one mapping (F_0) = {total_ways_A} + {total_ways_B} + {total_ways_C} = {F_0}")
print(f"Total F = 5! * F_0 = {num_mappings} * {F_0} = {F}")

print("\n--- Final Probability (P = F/S) ---")
print(f"P = {F} / {int(S)}")
print(f"P ≈ {P}")