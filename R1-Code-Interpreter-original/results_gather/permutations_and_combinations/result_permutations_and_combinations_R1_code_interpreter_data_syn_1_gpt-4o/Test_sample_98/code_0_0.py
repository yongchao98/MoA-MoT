from itertools import permutations

# Books to arrange in positions 1 to 6
books = ['B', 'C', 'D', 'E', 'F', 'G']

# Iterate over all permutations of the books
for perm in permutations(books):
    # Check the constraints
    if perm.index('G') != 2 and \
       perm.index('E') < perm.index('F') and \
       perm.index('G') < 5 and \
       perm.index('C') < perm.index('G') and \
       perm.index('E') > perm.index('B') and \
       perm.index('F') < perm.index('D'):
        # If all constraints are satisfied, print the arrangement
        arrangement = list(perm) + ['A']
        print(arrangement)
        break