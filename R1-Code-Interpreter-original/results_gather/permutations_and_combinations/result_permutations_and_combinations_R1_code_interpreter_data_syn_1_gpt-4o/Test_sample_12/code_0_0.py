from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if abs(arrangement.index('A') - arrangement.index('F')) != 1:
        return False
    if abs(arrangement.index('E') - arrangement.index('F')) != 1:
        return False
    if arrangement.index('B') == 1:
        return False
    if arrangement.index('F') == 0:
        return False
    if arrangement.index('F') < arrangement.index('C'):
        return False
    if arrangement.index('D') > arrangement.index('E'):
        return False
    if arrangement.index('C') == 1:
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break