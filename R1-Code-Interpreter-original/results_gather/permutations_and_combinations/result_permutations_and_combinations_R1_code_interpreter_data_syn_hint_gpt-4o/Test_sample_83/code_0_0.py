from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if abs(arrangement.index('B') - arrangement.index('C')) != 1:
        return False
    if arrangement.index('E') == 2:
        return False
    if abs(arrangement.index('C') - arrangement.index('E')) != 1:
        return False
    if arrangement.index('G') > arrangement.index('A'):
        return False
    if arrangement.index('G') == 6:
        return False
    if arrangement.index('E') > arrangement.index('F'):
        return False
    if arrangement.index('A') == 6:
        return False
    if arrangement.index('B') == 6:
        return False
    if arrangement.index('B') == 1:
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break