import math

# Total number of players
total_players = 12

# Number of players in a lineup
lineup_size = 5

# Calculate total combinations without restriction
total_combinations = math.comb(total_players, lineup_size)

# Calculate combinations including both Bob and Yogi
# If Bob and Yogi are in the lineup, we choose 3 more players from the remaining 10
combinations_with_both = math.comb(total_players - 2, lineup_size - 2)

# Valid combinations (total - those with both Bob and Yogi)
valid_combinations = total_combinations - combinations_with_both

print(valid_combinations)