from itertools import permutations

# Possible pieces for each round based on the analysis
round1_options = ['A', 'B', 'C']
round2_options = ['F', 'H']
round3_piece = 'G'
round4_options = ['B', 'C', 'D']

# Generate all possible combinations
possible_combinations = permutations(round1_options + round2_options + [round3_piece] + round4_options, 4)

# Check each combination to see if it satisfies the conditions
for combination in possible_combinations:
    if (combination[0] in round1_options and
        combination[1] in round2_options and
        combination[2] == round3_piece and
        combination[3] in round4_options):
        print(combination)
        break