from itertools import permutations

# Fixed position for C
books = ['A', 'B', 'D', 'E', 'F']
fixed_position = ['C']

# Generate all permutations of the remaining books
for perm in permutations(books):
    # Insert C at the beginning
    arrangement = fixed_position + list(perm)
    
    # Check constraints
    if (arrangement.index('B') == arrangement.index('E') + 1 or arrangement.index('B') == arrangement.index('E') - 1) and \
       (arrangement.index('A') == arrangement.index('F') + 1 or arrangement.index('A') == arrangement.index('F') - 1) and \
       arrangement.index('D') < arrangement.index('B') and \
       arrangement.index('A') < arrangement.index('D'):
        print(arrangement)
        break