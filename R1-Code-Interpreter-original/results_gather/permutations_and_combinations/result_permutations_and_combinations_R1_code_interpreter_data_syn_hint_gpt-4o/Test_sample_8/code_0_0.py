from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E']

# Iterate over all permutations of the books
for perm in permutations(books):
    # Check the constraints
    if perm.index('C') != 2 and perm.index('A') != 2:  # C and A cannot be in position 3
        if abs(perm.index('A') - perm.index('D')) == 1:  # D must be adjacent to A
            if perm.index('C') < 5:  # C must be to the left of F (position 6)
                if perm.index('A') < 5:  # A must be to the left of F (position 6)
                    if perm.index('B') < perm.index('E'):  # B must be to the left of E
                        # If all constraints are satisfied, print the arrangement
                        print(list(perm) + ['F'])
                        break