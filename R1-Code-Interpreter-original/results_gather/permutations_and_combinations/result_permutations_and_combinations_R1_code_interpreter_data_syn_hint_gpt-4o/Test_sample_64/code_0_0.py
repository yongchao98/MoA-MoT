from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Fixed position for G
fixed_positions = {6: 'G'}

# Function to check if a permutation satisfies all constraints
def is_valid(arrangement):
    # Check constraints
    if arrangement.index('D') not in [arrangement.index('F') - 1, arrangement.index('F') + 1]:
        return False
    if arrangement.index('A') > arrangement.index('D'):
        return False
    if arrangement.index('F') not in [arrangement.index('A') - 1, arrangement.index('A') + 1]:
        return False
    if arrangement.index('F') == 0:
        return False
    if arrangement.index('C') < arrangement.index('A'):
        return False
    if arrangement.index('E') == 3:
        return False
    if arrangement.index('A') > 5:
        return False
    return True

# Try all permutations of the books
for perm in permutations(books):
    # Insert G at the fixed position
    arrangement = list(perm[:6]) + ['G']
    if is_valid(arrangement):
        print(arrangement)
        break