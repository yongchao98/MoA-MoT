from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    # Constraint 1: Book A must be to the right of book E.
    if arrangement.index('A') <= arrangement.index('E'):
        return False
    # Constraint 2: Book G must be adjacent to book A.
    if abs(arrangement.index('G') - arrangement.index('A')) != 1:
        return False
    # Constraint 3: Book B must be placed in position 2.
    if arrangement[1] != 'B':
        return False
    # Constraint 4: Book E must be to the left of book C.
    if arrangement.index('E') >= arrangement.index('C'):
        return False
    # Constraint 5: Book C must be to the right of book E.
    if arrangement.index('C') <= arrangement.index('E'):
        return False
    # Constraint 6: Book E must be to the right of book B.
    if arrangement.index('E') <= arrangement.index('B'):
        return False
    # Constraint 7: Book D must be adjacent to book G.
    if abs(arrangement.index('D') - arrangement.index('G')) != 1:
        return False
    # Constraint 8: Book C cannot be placed in position 5.
    if arrangement[4] == 'C':
        return False
    # Constraint 9: Book A must be placed in position 4.
    if arrangement[3] != 'A':
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break