from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Fixed position for A
fixed_position = 3

# Generate all permutations of the other books
for perm in permutations([b for b in books if b != 'A']):
    # Insert A at the fixed position
    arrangement = list(perm)
    arrangement.insert(fixed_position - 1, 'A')
    
    # Check all constraints
    if (arrangement.index('F') in [arrangement.index('A') - 1, arrangement.index('A') + 1] and  # F adjacent to A
        arrangement.index('B') < arrangement.index('C') and  # B left of C
        arrangement.index('D') > arrangement.index('A') and  # D right of A
        arrangement.index('F') > arrangement.index('C') and  # F right of C
        arrangement.index('F') in [arrangement.index('E') - 1, arrangement.index('E') + 1]):  # F adjacent to E
        print(arrangement)
        break