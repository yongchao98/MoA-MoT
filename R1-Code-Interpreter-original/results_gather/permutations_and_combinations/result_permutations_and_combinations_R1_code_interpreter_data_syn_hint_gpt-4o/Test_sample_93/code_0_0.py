from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed positions
fixed_positions = {4: 'D'}

# Constraints
def is_valid(arrangement):
    # Check fixed positions
    if arrangement[4] != 'D':
        return False
    # Constraint 1: E cannot be in position 3
    if arrangement[2] == 'E':
        return False
    # Constraint 2: G cannot be in position 2
    if arrangement[1] == 'G':
        return False
    # Constraint 5: D must be to the left of G
    if arrangement.index('D') > arrangement.index('G'):
        return False
    # Constraint 6: B cannot be in position 6
    if arrangement[5] == 'B':
        return False
    # Constraint 7: D must be to the right of B
    if arrangement.index('D') < arrangement.index('B'):
        return False
    # Constraint 8: F cannot be in position 6
    if arrangement[5] == 'F':
        return False
    return True

# Generate all permutations and check constraints
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break