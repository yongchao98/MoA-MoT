import itertools

# Books to arrange
books = ['A', 'B', 'C', 'E', 'F']

# Fixed position for D
fixed_positions = {3: 'D'}

# Constraints
def is_valid(arrangement):
    # Check if C is not in position 1
    if arrangement[0] == 'C':
        return False
    # Check if E is not in position 6
    if arrangement[5] == 'E':
        return False
    # Check if A is not in position 5
    if arrangement[4] == 'A':
        return False
    # Check if B is not in position 3
    if arrangement[2] == 'B':
        return False
    # Check if F is to the left of B
    if arrangement.index('F') > arrangement.index('B'):
        return False
    # Check if A is adjacent to E
    if abs(arrangement.index('A') - arrangement.index('E')) != 1:
        return False
    # Check if A is to the left of C
    if arrangement.index('A') > arrangement.index('C'):
        return False
    return True

# Generate all permutations of the remaining books
for perm in itertools.permutations(books):
    # Insert D in the fixed position
    arrangement = list(perm)
    arrangement.insert(2, 'D')  # D is at position 3 (index 2)
    
    # Check if the arrangement is valid
    if is_valid(arrangement):
        print(arrangement)
        break