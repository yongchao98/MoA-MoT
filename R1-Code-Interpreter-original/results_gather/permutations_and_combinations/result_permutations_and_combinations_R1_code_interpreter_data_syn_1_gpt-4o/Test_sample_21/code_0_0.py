import itertools

# Books to arrange
books = ['A', 'C', 'D', 'E', 'F']

# Fixed position for B
fixed_position = 5

# Iterate over all permutations of the other books
for perm in itertools.permutations(books):
    # Create a list with B in the fixed position
    arrangement = list(perm[:fixed_position]) + ['B'] + list(perm[fixed_position:])
    
    # Check all constraints
    if (arrangement.index('F') < arrangement.index('C') and  # F is to the left of C
        arrangement.index('E') < arrangement.index('D') and  # E is to the left of D
        arrangement.index('A') < arrangement.index('E') and  # A is to the left of E
        arrangement.index('C') > arrangement.index('A') and  # C is to the right of A
        abs(arrangement.index('A') - arrangement.index('C')) == 1 and  # A is adjacent to C
        arrangement.index('F') < arrangement.index('A')):  # F is to the left of A
        print(arrangement)
        break