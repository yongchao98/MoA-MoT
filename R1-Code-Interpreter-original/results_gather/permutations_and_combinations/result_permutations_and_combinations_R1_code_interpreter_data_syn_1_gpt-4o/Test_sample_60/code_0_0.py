from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# E is fixed at position 6
fixed_positions = {5: 'E'}

def is_valid(arrangement):
    # Check all constraints
    if arrangement[4] != 'E':  # E must be at position 6 (index 5)
        return False
    if arrangement[6] == 'E':  # E cannot be at position 7
        return False
    if arrangement[4] == 'G':  # G cannot be at position 5
        return False
    if arrangement[5] == 'A':  # A cannot be at position 6
        return False
    if arrangement.index('G') > arrangement.index('D'):  # G must be to the left of D
        return False
    if arrangement.index('A') > arrangement.index('C'):  # A must be to the left of C
        return False
    if arrangement.index('A') > arrangement.index('G'):  # A must be to the left of G
        return False
    if arrangement.index('F') < arrangement.index('B'):  # F must be to the right of B
        return False
    if arrangement.index('F') < arrangement.index('D'):  # F must be to the right of D
        return False
    return True

# Generate all permutations of the books
for perm in permutations(books):
    # Apply fixed positions
    perm = list(perm)
    for pos, book in fixed_positions.items():
        if perm[pos] != book:
            break
    else:
        # Check if the current permutation is valid
        if is_valid(perm):
            print(perm)
            break