import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if abs(arrangement.index('C') - arrangement.index('E')) != 1:
        return False
    if arrangement.index('E') > arrangement.index('D'):
        return False
    if arrangement.index('A') > arrangement.index('D'):
        return False
    if arrangement.index('B') == 3:
        return False
    if arrangement.index('E') == 0:
        return False
    if abs(arrangement.index('C') - arrangement.index('B')) != 1:
        return False
    if arrangement.index('B') == 6:
        return False
    if arrangement.index('D') == 6:
        return False
    if arrangement.index('D') > arrangement.index('G'):
        return False
    return True

for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break