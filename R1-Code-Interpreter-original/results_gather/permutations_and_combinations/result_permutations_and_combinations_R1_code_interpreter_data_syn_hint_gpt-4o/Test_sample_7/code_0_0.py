from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('D') < arrangement.index('C'):
        return False
    if arrangement.index('F') > arrangement.index('A'):
        return False
    if arrangement.index('F') > arrangement.index('B'):
        return False
    if arrangement.index('C') == 4:
        return False
    if arrangement.index('A') == 3:
        return False
    if arrangement.index('E') == 3:
        return False
    if arrangement.index('C') == 2:
        return False
    if arrangement.index('F') > arrangement.index('D'):
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break