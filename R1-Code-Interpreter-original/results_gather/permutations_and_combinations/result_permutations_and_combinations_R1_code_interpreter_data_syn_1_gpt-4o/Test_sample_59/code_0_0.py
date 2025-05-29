from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement[3] != 'E':  # E must be in position 4
        return False
    if 'D' not in arrangement or 'F' not in arrangement:
        return False
    if abs(arrangement.index('D') - arrangement.index('F')) != 1:  # D and F must be adjacent
        return False
    if arrangement.index('F') < arrangement.index('D'):  # F must be to the right of D
        return False
    if arrangement.index('B') < arrangement.index('G'):  # B must be to the right of G
        return False
    if arrangement.index('B') == 1:  # B cannot be in position 2
        return False
    if arrangement.index('D') < arrangement.index('E'):  # D must be to the right of E
        return False
    if abs(arrangement.index('A') - arrangement.index('C')) != 1:  # A and C must be adjacent
        return False
    if arrangement.index('F') < arrangement.index('A'):  # F must be to the right of A
        return False
    return True

# Check all permutations
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break