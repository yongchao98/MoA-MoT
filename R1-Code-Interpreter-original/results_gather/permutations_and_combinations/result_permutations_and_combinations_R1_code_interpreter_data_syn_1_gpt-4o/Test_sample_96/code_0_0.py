import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    # Constraint 1: Book D cannot be in position 1
    if arrangement[0] == 'D':
        return False
    # Constraint 2: Book F cannot be in position 2
    if arrangement[1] == 'F':
        return False
    # Constraint 3: Book A must be to the left of book F
    if arrangement.index('A') > arrangement.index('F'):
        return False
    # Constraint 4: Book D cannot be in position 4
    if arrangement[3] == 'D':
        return False
    # Constraint 5: Book D cannot be in position 6
    if arrangement[5] == 'D':
        return False
    # Constraint 6: Book F cannot be in position 6
    if arrangement[5] == 'F':
        return False
    # Constraint 7: Book A cannot be in position 4
    if arrangement[3] == 'A':
        return False
    # Constraint 8: Book A must be to the left of book C
    if arrangement.index('A') > arrangement.index('C'):
        return False
    # Constraint 9: Book E must be to the left of book F
    if arrangement.index('E') > arrangement.index('F'):
        return False
    return True

# Find a valid arrangement
for arrangement in itertools.permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break