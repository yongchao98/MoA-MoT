from itertools import permutations

# Positions for books
positions = [1, 2, 3, 4, 5, 6, 7]

# Fixed positions
fixed_positions = {'A': 7, 'C': 5}

# Possible positions for B, E, F, G
possible_positions = [1, 2, 3, 4]

# Iterate over permutations of B, E, F
for perm in permutations(possible_positions, 3):
    B, E, F = perm
    # G must be to the right of B and F, so G is in position 6
    G = 6
    # Check constraints
    if E < 5 and B < G and F < G:
        # Create the arrangement
        arrangement = [''] * 7
        arrangement[B - 1] = 'B'
        arrangement[E - 1] = 'E'
        arrangement[F - 1] = 'F'
        arrangement[4] = 'C'  # C is in position 5
        arrangement[5] = 'G'  # G is in position 6
        arrangement[6] = 'A'  # A is in position 7
        print(arrangement)
        break