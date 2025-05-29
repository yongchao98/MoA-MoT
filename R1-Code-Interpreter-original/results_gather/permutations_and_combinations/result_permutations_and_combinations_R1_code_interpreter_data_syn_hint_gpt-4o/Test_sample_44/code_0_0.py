from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if abs(arrangement.index('A') - arrangement.index('D')) != 1:
        return False
    if arrangement.index('D') < arrangement.index('A'):
        return False
    if abs(arrangement.index('A') - arrangement.index('B')) != 1:
        return False
    if arrangement.index('A') > arrangement.index('C'):
        return False
    if arrangement.index('F') < arrangement.index('B'):
        return False
    if arrangement.index('E') == 2 or arrangement.index('E') == 5:
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break