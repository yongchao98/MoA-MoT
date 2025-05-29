import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    # Constraint 1: E must be to the left of F
    if arrangement.index('E') > arrangement.index('F'):
        return False
    # Constraint 2: G must be to the left of A
    if arrangement.index('G') > arrangement.index('A'):
        return False
    # Constraint 3: E must be adjacent to B
    if abs(arrangement.index('E') - arrangement.index('B')) != 1:
        return False
    # Constraint 4: D cannot be in position 7
    if arrangement.index('D') == 6:
        return False
    # Constraint 5: D must be adjacent to C
    if abs(arrangement.index('D') - arrangement.index('C')) != 1:
        return False
    # Constraint 6: G cannot be in position 1
    if arrangement.index('G') == 0:
        return False
    # Constraint 7: D must be to the right of C
    if arrangement.index('D') < arrangement.index('C'):
        return False
    # Constraint 8: B cannot be in position 7
    if arrangement.index('B') == 6:
        return False
    # Constraint 9: A must be to the right of G
    if arrangement.index('A') < arrangement.index('G'):
        return False
    return True

# Generate all permutations of the books
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break