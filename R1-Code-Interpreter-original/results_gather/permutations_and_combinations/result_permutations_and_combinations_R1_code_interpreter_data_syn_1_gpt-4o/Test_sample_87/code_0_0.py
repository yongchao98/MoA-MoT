from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for book F
fixed_position = 2

# Function to check if a permutation satisfies all constraints
def is_valid(arrangement):
    # Check if F is in position 2
    if arrangement[fixed_position - 1] != 'F':
        return False
    # Check constraints
    if arrangement.index('C') > arrangement.index('D'):
        return False
    if arrangement.index('A') > arrangement.index('C'):
        return False
    if arrangement.index('A') < arrangement.index('E'):
        return False
    if arrangement.index('C') < arrangement.index('B'):
        return False
    if arrangement.index('G') == 5 or arrangement.index('G') == 1:
        return False
    return True

# Try all permutations and find a valid one
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break