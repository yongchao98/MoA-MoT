from itertools import permutations

# Possible opponent pieces for each round
round1_options = {'A', 'B', 'C', 'D', 'E', 'F'}
round2_options = {'B', 'C', 'D', 'E', 'F', 'G'}
round3_options = {'A', 'B', 'C', 'D', 'E'}

# Find a valid combination of pieces
for perm in permutations(round1_options & round2_options & round3_options, 3):
    if perm[0] in round1_options and perm[1] in round2_options and perm[2] in round3_options:
        print(perm)
        break