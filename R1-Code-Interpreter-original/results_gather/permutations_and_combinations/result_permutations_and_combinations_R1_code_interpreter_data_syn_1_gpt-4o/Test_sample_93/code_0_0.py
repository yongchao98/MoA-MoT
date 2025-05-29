from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for D
fixed_positions = {4: 'D'}

# Constraints
def is_valid(arrangement):
    # Check all constraints
    if arrangement[2] == 'E':  # E cannot be in position 3
        return False
    if arrangement[1] == 'G':  # G cannot be in position 2
        return False
    if arrangement[4] != 'D':  # D must be in position 5
        return False
    if arrangement.index('D') > arrangement.index('G'):  # D must be to the left of G
        return False
    if arrangement.index('D') < arrangement.index('B'):  # D must be to the right of B
        return False
    if arrangement[5] == 'B':  # B cannot be in position 6
        return False
    if arrangement[5] == 'F':  # F cannot be in position 6
        return False
    return True

# Generate all permutations of the books
for perm in permutations(books):
    # Apply fixed positions
    if all(perm[pos] == book for pos, book in fixed_positions.items()):
        # Check if the arrangement is valid
        if is_valid(perm):
            print(list(perm))
            break