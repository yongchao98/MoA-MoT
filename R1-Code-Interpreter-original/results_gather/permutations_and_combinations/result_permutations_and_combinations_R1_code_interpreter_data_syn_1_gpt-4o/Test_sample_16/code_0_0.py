from itertools import permutations

# Fixed positions
positions = ['_', '_', 'A', 'D', '_', '_']

# Possible positions for B, F, C, E
possible_positions = [0, 1, 4, 5]

# Generate permutations of B, F, C, E
for perm in permutations(['B', 'F', 'C', 'E']):
    # Place B, F, C, E in the possible positions
    temp_positions = positions[:]
    for i, pos in enumerate(possible_positions):
        temp_positions[pos] = perm[i]
    
    # Check all constraints
    if (temp_positions[1] != 'F' and  # F cannot be in position 2
        temp_positions.index('C') > temp_positions.index('D') and  # C must be to the right of D
        abs(temp_positions.index('F') - temp_positions.index('B')) == 1 and  # F must be adjacent to B
        temp_positions.index('B') != 2):  # B cannot be in position 3
        print(temp_positions)
        break