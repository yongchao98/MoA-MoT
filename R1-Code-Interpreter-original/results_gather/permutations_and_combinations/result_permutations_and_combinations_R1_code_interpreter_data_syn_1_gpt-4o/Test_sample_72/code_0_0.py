import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    try:
        pos_A = arrangement.index('A')
        pos_B = arrangement.index('B')
        pos_C = arrangement.index('C')
        pos_D = arrangement.index('D')
        pos_E = arrangement.index('E')
        pos_F = arrangement.index('F')
        pos_G = arrangement.index('G')
        
        # Constraint 1: A must be adjacent to G
        if abs(pos_A - pos_G) != 1:
            return False
        
        # Constraint 2: C must be adjacent to F
        if abs(pos_C - pos_F) != 1:
            return False
        
        # Constraint 3: C must be to the left of G
        if pos_C >= pos_G:
            return False
        
        # Constraint 4: B must be to the left of G
        if pos_B >= pos_G:
            return False
        
        # Constraint 5: B must be adjacent to C
        if abs(pos_B - pos_C) != 1:
            return False
        
        # Constraint 6: F must be to the left of D
        if pos_F >= pos_D:
            return False
        
        # Constraint 7: B cannot be in position 1
        if pos_B == 0:
            return False
        
        # Constraint 8: G must be to the right of F
        if pos_G <= pos_F:
            return False
        
        return True
    except ValueError:
        return False

# Generate all permutations of the books
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break