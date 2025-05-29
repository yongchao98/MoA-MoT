from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('D') > arrangement.index('F'):
        return False
    if abs(arrangement.index('D') - arrangement.index('A')) != 1:
        return False
    if arrangement.index('F') > arrangement.index('C'):
        return False
    if arrangement.index('C') == 0:
        return False
    if arrangement.index('B') == 0:
        return False
    if arrangement.index('C') < arrangement.index('B'):
        return False
    if arrangement.index('A') > arrangement.index('C'):
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break