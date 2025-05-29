import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for G
fixed_position = 6

# Generate all permutations of the first 6 positions
for perm in itertools.permutations(books[:-1]):
    # Insert G at the fixed position
    arrangement = list(perm[:fixed_position]) + ['G'] + list(perm[fixed_position:])
    
    # Check all constraints
    if (arrangement.index('A') > arrangement.index('D') and
        arrangement.index('G') > arrangement.index('B') and
        arrangement.index('D') < arrangement.index('E') and
        arrangement.index('C') < arrangement.index('A') and
        abs(arrangement.index('D') - arrangement.index('E')) == 1 and
        arrangement.index('A') != 3 and
        arrangement.index('G') == 6 and
        arrangement.index('G') > arrangement.index('C')):
        print(arrangement)
        break