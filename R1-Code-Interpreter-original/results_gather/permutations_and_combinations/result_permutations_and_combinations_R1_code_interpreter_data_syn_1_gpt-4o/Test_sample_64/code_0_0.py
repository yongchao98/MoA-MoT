from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement[6] != 'G':  # G must be in position 7
        return False
    if 'D' not in [arrangement[i] for i in range(6) if arrangement[i] == 'F' and i < 6]:  # D must be adjacent to F
        return False
    if 'A' not in [arrangement[i] for i in range(6) if arrangement[i] == 'F' and i > 0]:  # F must be adjacent to A
        return False
    if arrangement.index('A') > arrangement.index('D'):  # A must be to the left of D
        return False
    if arrangement.index('A') > arrangement.index('C'):  # A must be to the left of C
        return False
    if arrangement.index('A') > arrangement.index('G'):  # A must be to the left of G
        return False
    if arrangement.index('F') == 0:  # F cannot be in position 1
        return False
    if arrangement.index('E') == 3:  # E cannot be in position 4
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break