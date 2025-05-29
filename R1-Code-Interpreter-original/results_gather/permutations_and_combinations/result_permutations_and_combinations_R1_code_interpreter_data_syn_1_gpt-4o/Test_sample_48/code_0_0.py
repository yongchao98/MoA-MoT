from itertools import permutations

# Define the books
books = ['A', 'B', 'C', 'D', 'E', 'F']

# Define the constraints
def is_valid(arrangement):
    # Constraint 1: Book D must be adjacent to book C
    if abs(arrangement.index('D') - arrangement.index('C')) != 1:
        return False
    # Constraint 2: Book A must be adjacent to book C
    if abs(arrangement.index('A') - arrangement.index('C')) != 1:
        return False
    # Constraint 3: Book A must be to the left of book B
    if arrangement.index('A') >= arrangement.index('B'):
        return False
    # Constraint 4: Book A must be to the left of book D
    if arrangement.index('A') >= arrangement.index('D'):
        return False
    # Constraint 5: Book A must be placed in position 2
    if arrangement[1] != 'A':
        return False
    # Constraint 6: Book D cannot be placed in position 1
    if arrangement[0] == 'D':
        return False
    # Constraint 7: Book B must be adjacent to book F
    if abs(arrangement.index('B') - arrangement.index('F')) != 1:
        return False
    # Constraint 8: Book A cannot be placed in position 1
    if arrangement[0] == 'A':
        return False
    return True

# Find a valid arrangement
for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break