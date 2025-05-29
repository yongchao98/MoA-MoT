import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']

# Generate all permutations of the books
for perm in itertools.permutations(books):
    # Check all constraints
    if perm.index('A') < perm.index('B') and \
       perm.index('D') != 0 and \
       perm.index('C') != 1 and \
       perm.index('F') != 1 and \
       abs(perm.index('D') - perm.index('A')) == 1 and \
       perm.index('D') != 5 and \
       perm.index('D') < perm.index('E') and \
       perm.index('C') != 4:
        print(list(perm))
        break