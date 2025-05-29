# Number of ways to choose apples (0 to 5)
ways_apples = 6

# Number of ways to choose oranges (0 to 10)
ways_oranges = 11

# Total combinations including the empty basket
total_combinations = ways_apples * ways_oranges

# Exclude the empty basket
total_baskets = total_combinations - 1

print(total_baskets)