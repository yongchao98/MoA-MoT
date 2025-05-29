from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement[5] != 'B':  # Book B must be in position 6
        return False
    if abs(arrangement.index('C') - arrangement.index('E')) != 1:  # C must be adjacent to E
        return False
    if arrangement.index('A') < arrangement.index('E'):  # A must be to the right of E
        return False
    if arrangement.index('E') < arrangement.index('F'):  # E must be to the right of F
        return False
    if arrangement.index('A') < arrangement.index('C'):  # A must be to the right of C
        return False
    if arrangement.index('G') < arrangement.index('E'):  # G must be to the right of E
        return False
    if arrangement.index('F') == 6:  # F cannot be in position 7
        return False
    if arrangement.index('C') > arrangement.index('G'):  # C must be to the left of G
        return False
    return True

# Try all permutations and find a valid one
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break