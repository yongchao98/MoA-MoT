from itertools import permutations

# Books to arrange
books = ['A', 'B', 'D', 'F', 'G']

# Fixed positions
fixed_positions = {2: 'E', 6: 'C'}

# Constraints
def is_valid(arrangement):
    # Place fixed books
    arrangement = list(arrangement)
    arrangement.insert(1, 'E')  # E is at position 2
    arrangement.insert(5, 'C')  # C is at position 6
    
    # Check constraints
    if arrangement.index('D') < arrangement.index('E'):
        return False
    if arrangement.index('G') < arrangement.index('A'):
        return False
    if arrangement.index('A') > arrangement.index('F'):
        return False
    if abs(arrangement.index('B') - arrangement.index('G')) != 1:
        return False
    if abs(arrangement.index('F') - arrangement.index('G')) != 1:
        return False
    if arrangement.index('F') == 6:
        return False
    return True

# Generate all permutations of the remaining books
for perm in permutations(books):
    if is_valid(perm):
        # Insert fixed books into the arrangement
        arrangement = list(perm)
        arrangement.insert(1, 'E')  # E is at position 2
        arrangement.insert(5, 'C')  # C is at position 6
        print(arrangement)
        break