import math

# Total number of ways to seat Mary and James
total_ways = math.comb(7, 2) * 2

# Number of ways to seat Mary and James next to each other
ways_next_to_each_other = 6 * 2

# Probability that they are not sitting next to each other
probability_not_next_to_each_other = 1 - (ways_next_to_each_other / total_ways)

print(probability_not_next_to_each_other)