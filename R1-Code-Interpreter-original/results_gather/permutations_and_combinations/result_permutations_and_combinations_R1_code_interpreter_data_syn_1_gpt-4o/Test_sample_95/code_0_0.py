import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    # Constraint 1: Book C cannot be placed in position 3
    if arrangement[2] == 'C':
        return False
    # Constraint 2: Book A cannot be placed in position 1
    if arrangement[0] == 'A':
        return False
    # Constraint 3: Book C cannot be placed in position 5
    if arrangement[4] == 'C':
        return False
    # Constraint 4: Book E must be placed in position 5
    if arrangement[4] != 'E':
        return False
    # Constraint 5: Book E must be to the left of book F
    if arrangement.index('E') > arrangement.index('F'):
        return False
    # Constraint 6: Book A cannot be placed in position 6
    if arrangement[5] == 'A':
        return False
    # Constraint 7: Book D must be adjacent to book B
    if abs(arrangement.index('D') - arrangement.index('B')) != 1:
        return False
    # Constraint 8: Book B must be to the left of book C
    if arrangement.index('B') > arrangement.index('C'):
        return False
    # Constraint 9: Book D must be to the left of book C
    if arrangement.index('D') > arrangement.index('C'):
        return False
    return True

# Generate all permutations of the books
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break