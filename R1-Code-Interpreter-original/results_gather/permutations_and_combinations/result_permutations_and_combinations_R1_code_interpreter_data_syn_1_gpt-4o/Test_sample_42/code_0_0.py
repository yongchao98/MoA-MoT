import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if abs(arrangement.index('B') - arrangement.index('F')) != 1:
        return False
    if abs(arrangement.index('A') - arrangement.index('B')) != 1:
        return False
    if arrangement.index('F') == 2:
        return False
    if arrangement.index('D') < arrangement.index('E'):
        return False
    if arrangement.index('C') == 4:
        return False
    if abs(arrangement.index('F') - arrangement.index('C')) != 1:
        return False
    if arrangement.index('B') < arrangement.index('E'):
        return False
    return True

for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break