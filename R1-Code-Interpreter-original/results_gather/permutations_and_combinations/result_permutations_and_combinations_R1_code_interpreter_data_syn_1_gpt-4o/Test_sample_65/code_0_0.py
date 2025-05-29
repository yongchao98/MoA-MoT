import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('A') < arrangement.index('B'):
        return False
    if arrangement.index('F') < arrangement.index('G'):
        return False
    if abs(arrangement.index('D') - arrangement.index('E')) != 1:
        return False
    if arrangement.index('F') == 3:
        return False
    if arrangement.index('E') > arrangement.index('A'):
        return False
    if arrangement.index('C') == 3:
        return False
    if arrangement.index('B') > arrangement.index('G'):
        return False
    if arrangement.index('C') < arrangement.index('F'):
        return False
    return True

for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break