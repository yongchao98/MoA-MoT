from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    # Constraint 3: Book B must be placed in position 4
    if arrangement[3] != 'B':
        return False
    # Constraint 6: Book B cannot be placed in position 6
    if arrangement[5] == 'B':
        return False
    # Constraint 8: Book C cannot be placed in position 7
    if arrangement[6] == 'C':
        return False
    # Constraint 1: Book D must be to the left of book F
    if arrangement.index('D') > arrangement.index('F'):
        return False
    # Constraint 2: Book F must be to the right of book B
    if arrangement.index('F') < arrangement.index('B'):
        return False
    # Constraint 4: Book D cannot be placed in position 5
    if arrangement[4] == 'D':
        return False
    # Constraint 5: Book A must be to the left of book E
    if arrangement.index('A') > arrangement.index('E'):
        return False
    # Constraint 7: Book E must be to the left of book D
    if arrangement.index('E') > arrangement.index('D'):
        return False
    # Constraint 9: Book B must be to the right of book A
    if arrangement.index('B') < arrangement.index('A'):
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break