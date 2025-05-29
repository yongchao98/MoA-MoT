from itertools import product

# Possible pieces for each round based on the results
round1_possibilities = ['E', 'F', 'G', 'H']
round2_possibilities = ['B', 'C', 'D', 'E', 'F', 'G', 'H']
round3_possibilities = ['A', 'B', 'C', 'D']
round4_possibilities = ['A', 'B', 'C', 'D']

# Generate all possible combinations
all_combinations = product(round1_possibilities, round2_possibilities, round3_possibilities, round4_possibilities)

# Check for valid combinations
valid_combinations = []
for combo in all_combinations:
    if len(set(combo)) == 4:  # Ensure all pieces are unique
        valid_combinations.append(combo)

# Print the valid combinations
print(valid_combinations)