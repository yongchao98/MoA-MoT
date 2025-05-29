# Possible pieces for each round based on the analysis
round_1_possibilities = {'A', 'B', 'C', 'D', 'E'}
round_2_possibilities = {'H'}
round_3_possibilities = {'E', 'F', 'G'}
round_4_possibilities = {'B', 'C', 'D', 'E', 'F', 'G'}

# Since H is used in round 2, remove it from other possibilities
round_1_possibilities.discard('H')
round_3_possibilities.discard('H')
round_4_possibilities.discard('H')

# Deduce the exact pieces
# Round 3 must be E, F, or G, and it cannot be E if E is used in round 1
# Round 4 must be B, C, D, E, F, or G, and it cannot be E if E is used in round 1

# Check combinations
for r1 in round_1_possibilities:
    for r3 in round_3_possibilities:
        for r4 in round_4_possibilities:
            if r1 != r3 and r1 != r4 and r3 != r4:
                print([r1, 'H', r3, r4])