from math import comb

# Calculate the number of favorable outcomes
favorable_outcomes = (4 * 4) + (4 * 4) + (4 * 4) + (4 * 4) + (4 * 3)

# Calculate the total number of ways to choose 2 cards from 52
total_outcomes = comb(52, 2)

# Calculate the probability
probability = favorable_outcomes / total_outcomes

print(probability)