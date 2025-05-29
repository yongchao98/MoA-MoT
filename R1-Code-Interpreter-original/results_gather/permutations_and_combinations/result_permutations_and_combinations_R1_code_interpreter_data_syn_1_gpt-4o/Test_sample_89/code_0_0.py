import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Fixed position for C
fixed_position_C = 4

# Generate all permutations of the books
for perm in itertools.permutations(books):
    # Check if C is in the correct position
    if perm[fixed_position_C - 1] != 'C':
        continue
    
    # Extract positions of each book
    pos_B = perm.index('B') + 1
    pos_E = perm.index('E') + 1
    pos_D = perm.index('D') + 1
    pos_F = perm.index('F') + 1
    pos_A = perm.index('A') + 1
    
    # Check all constraints
    if (pos_B < pos_C and
        abs(pos_B - pos_E) == 1 and
        pos_D > pos_E and
        pos_F != 5 and
        pos_B != 2 and
        pos_F > pos_C and
        pos_A != 6):
        # If all constraints are satisfied, print the valid arrangement
        print(list(perm))
        break