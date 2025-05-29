import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Fixed position for F
fixed_position = 3

# Generate all permutations of the other books
other_books = ['A', 'B', 'C', 'D', 'E']
permutations = itertools.permutations(other_books)

# Check each permutation against the constraints
for perm in permutations:
    # Insert F at the fixed position
    arrangement = list(perm)
    arrangement.insert(fixed_position - 1, 'F')
    
    # Check constraints
    if (arrangement.index('B') > arrangement.index('E') and  # B is to the right of E
        arrangement.index('A') > arrangement.index('D') and  # A is to the right of D
        arrangement.index('E') != 0 and                      # E is not in position 1
        arrangement.index('C') != 1 and                      # C is not in position 2
        arrangement.index('B') != 1):                        # B is not in position 2
        print(arrangement)
        break