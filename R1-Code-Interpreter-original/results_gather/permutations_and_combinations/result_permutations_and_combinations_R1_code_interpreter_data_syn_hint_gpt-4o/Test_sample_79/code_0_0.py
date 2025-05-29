from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('A') < arrangement.index('B'):
        return False
    if arrangement.index('F') < arrangement.index('A'):
        return False
    if arrangement.index('F') == 2:
        return False
    if abs(arrangement.index('A') - arrangement.index('D')) != 1:
        return False
    if arrangement.index('C') == 4:
        return False
    if arrangement.index('B') > arrangement.index('F'):
        return False
    if abs(arrangement.index('A') - arrangement.index('C')) != 1:
        return False
    if arrangement.index('E') == 3:
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break