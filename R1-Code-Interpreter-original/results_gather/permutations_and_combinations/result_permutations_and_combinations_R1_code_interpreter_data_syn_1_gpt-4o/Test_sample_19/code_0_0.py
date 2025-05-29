from itertools import permutations

# Initial positions
positions = ['B', 'A', None, None, None, None]

# Possible books for remaining positions
remaining_books = ['C', 'D', 'E', 'F']

# Check all permutations of the remaining books
for perm in permutations(remaining_books):
    # Place the books in the remaining positions
    temp_positions = positions[:]
    temp_positions[2:] = perm
    
    # Check constraints
    if temp_positions.index('C') > temp_positions.index('B') and \
       temp_positions.index('F') > temp_positions.index('B') and \
       temp_positions.index('E') < temp_positions.index('F') and \
       temp_positions.index('D') != 2:
        print(temp_positions)
        break