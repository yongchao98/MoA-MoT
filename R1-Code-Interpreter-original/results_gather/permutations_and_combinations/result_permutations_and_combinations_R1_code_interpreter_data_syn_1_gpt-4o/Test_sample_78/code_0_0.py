import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for B
fixed_positions = {1: 'B'}

# Constraints
def is_valid(arrangement):
    # Check fixed position for B
    if arrangement[1] != 'B':
        return False
    # Check constraints
    if arrangement[3] == 'A':  # A cannot be in position 4
        return False
    if arrangement[2] == 'G':  # G cannot be in position 3
        return False
    if arrangement[4] == 'F':  # F cannot be in position 5
        return False
    if arrangement[2] == 'C':  # C cannot be in position 3
        return False
    if 'A' not in arrangement or 'F' not in arrangement:
        return False
    if abs(arrangement.index('A') - arrangement.index('F')) != 1:  # A must be adjacent to F
        return False
    if 'A' not in arrangement or 'E' not in arrangement:
        return False
    if abs(arrangement.index('A') - arrangement.index('E')) != 1:  # E must be adjacent to A
        return False
    if arrangement[3] == 'E':  # E cannot be in position 4
        return False
    return True

# Generate all permutations of the books except B
for perm in itertools.permutations([book for book in books if book != 'B']):
    # Insert B at position 2
    arrangement = list(perm)
    arrangement.insert(1, 'B')
    # Check if the arrangement is valid
    if is_valid(arrangement):
        print(arrangement)
        break