from itertools import permutations

# Positions 1, 2, 4, 5 are available for A, B, C, D
positions = [1, 2, 4, 5]
books = ['A', 'B', 'C', 'D']

# Constraint: B must be to the right of D, C must be to the right of B
for perm in permutations(positions):
    d_pos, b_pos, c_pos, a_pos = perm
    if d_pos < b_pos < c_pos:
        # Construct the full arrangement
        arrangement = [''] * 6
        arrangement[0] = books[perm.index(1)]
        arrangement[1] = books[perm.index(2)]
        arrangement[2] = 'F'  # F is in position 3
        arrangement[3] = books[perm.index(4)]
        arrangement[4] = books[perm.index(5)]
        arrangement[5] = 'E'  # E is in position 6
        print(arrangement)
        break