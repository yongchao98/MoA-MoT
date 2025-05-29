import math

# Total number of arrangements without considering symmetries
total_arrangements = math.factorial(10)

# Number of symmetries
num_symmetries = 10

# Number of distinct arrangements using Burnside's Lemma
distinct_arrangements = total_arrangements // num_symmetries

print(distinct_arrangements)