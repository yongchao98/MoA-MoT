import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('A') < arrangement.index('D'):
        return False
    if arrangement.index('F') == 0:
        return False
    if arrangement.index('G') == 4:
        return False
    if abs(arrangement.index('B') - arrangement.index('F')) != 1:
        return False
    if arrangement.index('B') == 5:
        return False
    if arrangement.index('D') > arrangement.index('B'):
        return False
    if arrangement.index('G') < arrangement.index('D'):
        return False
    if arrangement.index('E') < arrangement.index('G'):
        return False
    if abs(arrangement.index('F') - arrangement.index('A')) != 1:
        return False
    return True

# Generate all permutations of the books
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break