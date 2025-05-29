import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']

# Constraints
def is_valid(arrangement):
    # Constraint 1: A must be adjacent to F
    if abs(arrangement.index('A') - arrangement.index('F')) != 1:
        return False
    # Constraint 2: E must be adjacent to F
    if abs(arrangement.index('E') - arrangement.index('F')) != 1:
        return False
    # Constraint 3: B cannot be in position 2
    if arrangement[1] == 'B':
        return False
    # Constraint 4: F cannot be in position 1
    if arrangement[0] == 'F':
        return False
    # Constraint 5: F must be to the right of C
    if arrangement.index('F') < arrangement.index('C'):
        return False
    # Constraint 6: D must be to the left of E
    if arrangement.index('D') > arrangement.index('E'):
        return False
    # Constraint 7: C cannot be in position 2
    if arrangement[1] == 'C':
        return False
    return True

# Generate all permutations and check for validity
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break