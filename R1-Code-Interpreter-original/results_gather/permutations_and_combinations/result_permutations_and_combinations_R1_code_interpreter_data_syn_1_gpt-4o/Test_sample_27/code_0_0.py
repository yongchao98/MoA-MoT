from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Fixed position for B
fixed_position = 2

# Generate all permutations of the remaining books
for perm in permutations([book for book in books if book != 'B']):
    # Insert B at the fixed position
    arrangement = list(perm)
    arrangement.insert(fixed_position - 1, 'B')
    
    # Check all constraints
    if (arrangement.index('C') < arrangement.index('A') and  # C is to the left of A
        arrangement.index('E') != 5 and  # E is not in position 6
        arrangement.index('A') != 4 and  # A is not in position 5
        abs(arrangement.index('E') - arrangement.index('F')) == 1 and  # E is adjacent to F
        arrangement.index('B') < arrangement.index('D') and  # B is to the left of D
        abs(arrangement.index('F') - arrangement.index('B')) == 1):  # F is adjacent to B
        print(arrangement)
        break