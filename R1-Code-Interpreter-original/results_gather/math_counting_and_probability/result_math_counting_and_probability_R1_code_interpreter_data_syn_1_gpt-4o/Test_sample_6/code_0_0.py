import math

# Total possible outcomes
total_outcomes = 6 ** 6

# Favorable outcomes (permutations of 6 distinct numbers)
favorable_outcomes = math.factorial(6)

# Probability
probability = favorable_outcomes / total_outcomes

print(probability)