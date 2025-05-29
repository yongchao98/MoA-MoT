from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Constraints
def is_valid(arrangement):
    # Constraint 1: Book D must be in position 6
    if arrangement[5] != 'D':
        return False
    # Constraint 2: Book D must be to the right of book E
    if arrangement.index('D') < arrangement.index('E'):
        return False
    # Constraint 3: Book C cannot be in position 2
    if arrangement[1] == 'C':
        return False
    # Constraint 4: Book C must be adjacent to book B
    if abs(arrangement.index('C') - arrangement.index('B')) != 1:
        return False
    # Constraint 5: Book C must be to the left of book D
    if arrangement.index('C') > arrangement.index('D'):
        return False
    # Constraint 6: Book F must be adjacent to book B
    if abs(arrangement.index('F') - arrangement.index('B')) != 1:
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break