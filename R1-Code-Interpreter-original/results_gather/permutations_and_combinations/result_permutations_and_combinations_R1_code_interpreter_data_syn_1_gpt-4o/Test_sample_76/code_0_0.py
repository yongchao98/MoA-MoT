from itertools import permutations

# Fixed position for F
fixed_positions = {3: 'F'}

# Constraints
def is_valid(arrangement):
    # Check fixed position for F
    if arrangement[3] != 'F':
        return False
    # Check B is to the right of F
    if arrangement.index('B') <= 3:
        return False
    # Check A is not in position 7
    if arrangement[6] == 'A':
        return False
    # Check E is not in position 2 or 6
    if arrangement[1] == 'E' or arrangement[5] == 'E':
        return False
    # Check A is not in position 4
    if arrangement[3] == 'A':
        return False
    # Check A is to the left of C
    if arrangement.index('A') >= arrangement.index('C'):
        return False
    # Check A is to the right of F
    if arrangement.index('A') <= 3:
        return False
    # Check A is to the right of D
    if arrangement.index('A') <= arrangement.index('D'):
        return False
    return True

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Generate all permutations of books
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break