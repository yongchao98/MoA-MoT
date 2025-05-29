import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('D') == 6:
        return False
    if arrangement.index('A') > arrangement.index('F'):
        return False
    if arrangement.index('G') > arrangement.index('F'):
        return False
    if abs(arrangement.index('G') - arrangement.index('B')) != 1:
        return False
    if arrangement.index('C') > arrangement.index('G'):
        return False
    if arrangement.index('G') < arrangement.index('D'):
        return False
    if arrangement.index('E') < arrangement.index('C'):
        return False
    if arrangement.index('F') < arrangement.index('C'):
        return False
    if arrangement.index('G') < arrangement.index('E'):
        return False
    return True

for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break