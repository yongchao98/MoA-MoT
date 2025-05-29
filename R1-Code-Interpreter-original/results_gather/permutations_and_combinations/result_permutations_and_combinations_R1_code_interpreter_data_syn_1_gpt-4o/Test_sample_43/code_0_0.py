from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: Book D cannot be placed in position 6
    if arrangement[5] == 'D':
        return False
    # Constraint 2: Book A cannot be placed in position 2
    if arrangement[1] == 'A':
        return False
    # Constraint 3: Book E must be adjacent to book C
    if abs(arrangement.index('E') - arrangement.index('C')) != 1:
        return False
    # Constraint 4: Book D must be to the left of book A
    if arrangement.index('D') > arrangement.index('A'):
        return False
    # Constraint 5: Book C must be to the left of book A
    if arrangement.index('C') > arrangement.index('A'):
        return False
    # Constraint 6: Book A must be placed in position 6
    if arrangement[5] != 'A':
        return False
    # Constraint 7: Book C must be adjacent to book B
    if abs(arrangement.index('C') - arrangement.index('B')) != 1:
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break