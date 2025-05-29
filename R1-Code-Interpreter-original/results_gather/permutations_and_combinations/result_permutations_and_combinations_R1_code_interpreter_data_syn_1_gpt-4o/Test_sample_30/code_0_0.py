import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Constraints
def is_valid(arrangement):
    # Constraint 1: A must be to the left of F
    if arrangement.index('A') >= arrangement.index('F'):
        return False
    # Constraint 2: F must be to the left of D
    if arrangement.index('F') >= arrangement.index('D'):
        return False
    # Constraint 3: C must be to the left of B
    if arrangement.index('C') >= arrangement.index('B'):
        return False
    # Constraint 4: B must be to the right of A
    if arrangement.index('B') <= arrangement.index('A'):
        return False
    # Constraint 5: E must be to the right of B
    if arrangement.index('E') <= arrangement.index('B'):
        return False
    # Constraint 6: C cannot be in position 1
    if arrangement.index('C') == 0:
        return False
    return True

# Find a valid arrangement
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break