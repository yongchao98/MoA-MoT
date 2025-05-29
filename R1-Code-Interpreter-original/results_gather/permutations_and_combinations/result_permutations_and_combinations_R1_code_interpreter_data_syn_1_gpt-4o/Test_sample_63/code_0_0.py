import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for F
fixed_positions = {3: 'F'}

# Constraints
def is_valid(arrangement):
    # Constraint 1: G cannot be in position 7
    if arrangement[6] == 'G':
        return False
    # Constraint 4: E cannot be in position 2
    if arrangement[1] == 'E':
        return False
    # Constraint 6: D cannot be in position 7
    if arrangement[6] == 'D':
        return False
    # Constraint 2: F must be to the left of B
    if arrangement.index('F') > arrangement.index('B'):
        return False
    # Constraint 3: B must be to the right of D
    if arrangement.index('B') < arrangement.index('D'):
        return False
    # Constraint 5: B must be to the right of G
    if arrangement.index('B') < arrangement.index('G'):
        return False
    # Constraint 7: D must be to the left of G
    if arrangement.index('D') > arrangement.index('G'):
        return False
    return True

# Generate all permutations of the books
for perm in itertools.permutations(books):
    # Apply fixed positions
    if perm[2] == 'F' and is_valid(perm):
        print(list(perm))
        break