from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed positions
fixed_positions = {1: 'C', 4: 'D'}

# Function to check if a permutation satisfies all constraints
def is_valid(arrangement):
    # Check fixed positions
    if arrangement[0] != 'C' or arrangement[3] != 'D':
        return False
    # Check B is not in position 1
    if arrangement[0] == 'B':
        return False
    # Check F is to the right of D
    if arrangement.index('F') < arrangement.index('D'):
        return False
    # Check F is to the left of E
    if arrangement.index('F') > arrangement.index('E'):
        return False
    # Check G is adjacent to C
    if abs(arrangement.index('G') - arrangement.index('C')) != 1:
        return False
    return True

# Generate all permutations of the books
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break