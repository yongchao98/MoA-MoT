from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Fixed position for C
fixed_position = 2

# Generate all permutations of the remaining books
for perm in permutations([b for b in books if b != 'C']):
    # Insert C at the fixed position
    arrangement = list(perm)
    arrangement.insert(fixed_position - 1, 'C')
    
    # Check all constraints
    if (arrangement.index('C') < arrangement.index('F') and  # C is to the left of F
        abs(arrangement.index('B') - arrangement.index('F')) == 1 and  # B is adjacent to F
        abs(arrangement.index('F') - arrangement.index('E')) == 1 and  # F is adjacent to E
        arrangement.index('C') < arrangement.index('E') and  # E is to the right of C
        arrangement.index('A') < arrangement.index('F') and  # A is to the left of F
        abs(arrangement.index('D') - arrangement.index('C')) == 1):  # D is adjacent to C
        print(arrangement)
        break