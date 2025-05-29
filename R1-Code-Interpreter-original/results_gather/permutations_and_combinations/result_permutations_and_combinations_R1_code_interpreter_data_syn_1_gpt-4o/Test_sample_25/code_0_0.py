from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Constraints
def is_valid(arrangement):
    # Constraint 1: B must be to the right of C
    if arrangement.index('B') <= arrangement.index('C'):
        return False
    # Constraint 2: D must be to the left of E
    if arrangement.index('D') >= arrangement.index('E'):
        return False
    # Constraint 3: A cannot be in position 2
    if arrangement.index('A') == 1:
        return False
    # Constraint 4: C cannot be in position 1
    if arrangement.index('C') == 0:
        return False
    # Constraint 5: B must be to the left of D
    if arrangement.index('B') >= arrangement.index('D'):
        return False
    # Constraint 6: C must be to the left of D
    if arrangement.index('C') >= arrangement.index('D'):
        return False
    # Constraint 7: F must be adjacent to A
    if abs(arrangement.index('F') - arrangement.index('A')) != 1:
        return False
    return True

# Find a valid arrangement
for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break