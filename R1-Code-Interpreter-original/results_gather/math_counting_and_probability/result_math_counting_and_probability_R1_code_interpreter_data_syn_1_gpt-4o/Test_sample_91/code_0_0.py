import math

# Total players to choose from (excluding the triplets)
n = 11

# Number of players to choose
k = 3

# Calculate the number of combinations
combinations = math.comb(n, k)

print(combinations)