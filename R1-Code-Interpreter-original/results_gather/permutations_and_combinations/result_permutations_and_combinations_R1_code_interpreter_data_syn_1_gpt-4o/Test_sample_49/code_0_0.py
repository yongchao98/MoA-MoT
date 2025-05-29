from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('C') == 5:  # C cannot be in position 6
        return False
    if arrangement.index('B') < arrangement.index('F'):  # B must be to the right of F
        return False
    if arrangement.index('A') > arrangement.index('D'):  # A must be to the left of D
        return False
    if arrangement.index('F') < arrangement.index('E'):  # F must be to the right of E
        return False
    if arrangement.index('A') < arrangement.index('F'):  # A must be to the right of F
        return False
    if arrangement.index('B') < arrangement.index('E'):  # B must be to the right of E
        return False
    if arrangement.index('A') == 4:  # A cannot be in position 5
        return False
    if arrangement.index('B') == 1:  # B cannot be in position 2
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break