from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for B
fixed_positions = {1: 'B'}

# Constraints
def is_valid(arrangement):
    # Check fixed position for B
    if arrangement[1] != 'B':
        return False
    # Check position constraints
    if arrangement[3] == 'A' or arrangement[3] == 'E':
        return False
    if arrangement[4] == 'F':
        return False
    if arrangement[2] == 'G' or arrangement[2] == 'C':
        return False
    # Check adjacency constraints
    a_index = arrangement.index('A')
    f_index = arrangement.index('F')
    e_index = arrangement.index('E')
    if abs(a_index - f_index) != 1:
        return False
    if abs(a_index - e_index) != 1:
        return False
    return True

# Generate all permutations of the books
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break