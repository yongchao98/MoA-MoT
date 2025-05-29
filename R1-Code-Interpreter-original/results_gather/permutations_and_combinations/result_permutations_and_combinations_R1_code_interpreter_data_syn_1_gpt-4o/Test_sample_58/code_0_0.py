from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed positions
fixed_positions = {'B': 3, 'G': 2}

# Constraints
def is_valid(arrangement):
    # Check fixed positions
    if arrangement[1] != 'G' or arrangement[2] != 'B':
        return False
    # F must be adjacent to D
    if not (abs(arrangement.index('F') - arrangement.index('D')) == 1):
        return False
    # C must be adjacent to F
    if not (abs(arrangement.index('C') - arrangement.index('F')) == 1):
        return False
    # F cannot be in position 1
    if arrangement[0] == 'F':
        return False
    # E must be to the left of C
    if arrangement.index('E') > arrangement.index('C'):
        return False
    # A cannot be in position 4
    if arrangement[3] == 'A':
        return False
    # A must be to the left of D
    if arrangement.index('A') > arrangement.index('D'):
        return False
    return True

# Generate all permutations and check constraints
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break