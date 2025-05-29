from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('A') == 1:
        return False
    if arrangement.index('E') > arrangement.index('F'):
        return False
    if abs(arrangement.index('D') - arrangement.index('F')) != 1:
        return False
    if arrangement.index('F') < arrangement.index('C'):
        return False
    if arrangement.index('E') < arrangement.index('A'):
        return False
    if arrangement.index('G') > arrangement.index('B'):
        return False
    if arrangement.index('E') < arrangement.index('B'):
        return False
    if arrangement.index('F') == 5:
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break