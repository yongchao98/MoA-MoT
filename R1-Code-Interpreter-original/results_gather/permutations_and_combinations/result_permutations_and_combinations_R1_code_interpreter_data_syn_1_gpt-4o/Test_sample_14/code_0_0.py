from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'E', 'F']

# Fixed position for D
fixed_position = ['D']

# Iterate over all permutations of the other books
for perm in permutations(books):
    # Create a full arrangement with D in position 6
    arrangement = list(perm) + fixed_position
    
    # Check all constraints
    if (arrangement.index('C') != 1 and  # C cannot be in position 2
        arrangement.index('C') < arrangement.index('D') and  # C must be to the left of D
        abs(arrangement.index('C') - arrangement.index('B')) == 1 and  # C must be adjacent to B
        abs(arrangement.index('F') - arrangement.index('B')) == 1 and  # F must be adjacent to B
        arrangement.index('E') < arrangement.index('D')):  # E must be to the left of D
        print(arrangement)
        break