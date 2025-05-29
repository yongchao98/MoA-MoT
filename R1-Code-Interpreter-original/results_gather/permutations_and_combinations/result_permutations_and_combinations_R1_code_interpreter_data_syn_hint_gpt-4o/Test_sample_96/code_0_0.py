from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    pos_A = arrangement.index('A')
    pos_B = arrangement.index('B')
    pos_C = arrangement.index('C')
    pos_D = arrangement.index('D')
    pos_E = arrangement.index('E')
    pos_F = arrangement.index('F')
    pos_G = arrangement.index('G')
    
    # Check all constraints
    if pos_D == 0:  # D cannot be in position 1
        return False
    if pos_F == 1:  # F cannot be in position 2
        return False
    if pos_A >= pos_F:  # A must be to the left of F
        return False
    if pos_D == 3:  # D cannot be in position 4
        return False
    if pos_D == 5:  # D cannot be in position 6
        return False
    if pos_F == 5:  # F cannot be in position 6
        return False
    if pos_A == 3:  # A cannot be in position 4
        return False
    if pos_A >= pos_C:  # A must be to the left of C
        return False
    if pos_E >= pos_F:  # E must be to the left of F
        return False
    
    return True

# Find a valid arrangement
for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break