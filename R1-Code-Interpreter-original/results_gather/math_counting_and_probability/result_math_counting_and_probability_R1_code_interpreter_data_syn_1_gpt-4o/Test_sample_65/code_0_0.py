from itertools import permutations

# Define the possible combinations that sum to 16
combinations = [(6, 6, 4), (6, 4, 6), (4, 6, 6), (5, 5, 6), (5, 6, 5), (6, 5, 5)]

# Use a set to avoid counting duplicate permutations
unique_permutations = set()

# Calculate all unique permutations of each combination
for combo in combinations:
    unique_permutations.update(permutations(combo))

# Total number of favorable outcomes
favorable_outcomes = len(unique_permutations)

# Total possible outcomes when rolling three dice
total_outcomes = 6 ** 3

# Calculate the probability
probability = favorable_outcomes / total_outcomes

print(probability)