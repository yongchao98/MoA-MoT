from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Constraint 1: Book B cannot be placed in position 5.
    if arrangement[4] == 'B':
        return False
    # Constraint 2: Book E cannot be placed in position 7.
    if arrangement[6] == 'E':
        return False
    # Constraint 3: Book A must be adjacent to book B.
    if abs(arrangement.index('A') - arrangement.index('B')) != 1:
        return False
    # Constraint 4: Book G must be to the left of book D.
    if arrangement.index('G') > arrangement.index('D'):
        return False
    # Constraint 5: Book E must be to the right of book G.
    if arrangement.index('E') < arrangement.index('G'):
        return False
    # Constraint 6: Book G must be to the left of book C.
    if arrangement.index('G') > arrangement.index('C'):
        return False
    # Constraint 7: Book C cannot be placed in position 3.
    if arrangement[2] == 'C':
        return False
    # Constraint 8: Book D must be placed in position 5.
    if arrangement[4] != 'D':
        return False
    # Constraint 9: Book B must be to the right of book E.
    if arrangement.index('B') < arrangement.index('E'):
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break