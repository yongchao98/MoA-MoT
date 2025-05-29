import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Constraints
def is_valid(arrangement):
    # Constraint 1: B must be adjacent to E
    if abs(arrangement.index('B') - arrangement.index('E')) != 1:
        return False
    # Constraint 2: D cannot be in position 6
    if arrangement.index('D') == 5:
        return False
    # Constraint 3: A must be to the right of C
    if arrangement.index('A') <= arrangement.index('C'):
        return False
    # Constraint 4: D must be adjacent to A
    if abs(arrangement.index('D') - arrangement.index('A')) != 1:
        return False
    # Constraint 5: C must be to the right of F
    if arrangement.index('C') <= arrangement.index('F'):
        return False
    # Constraint 6: F must be adjacent to C
    if abs(arrangement.index('F') - arrangement.index('C')) != 1:
        return False
    # Constraint 7: A cannot be in position 1
    if arrangement.index('A') == 0:
        return False
    return True

# Generate all permutations of the books
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break