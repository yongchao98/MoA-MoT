from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: Book C cannot be placed in position 2
    if arrangement[1] == 'C':
        return False
    # Constraint 2: Book D must be to the left of book F
    if arrangement.index('D') > arrangement.index('F'):
        return False
    # Constraint 3: Book E must be adjacent to book B
    if abs(arrangement.index('E') - arrangement.index('B')) != 1:
        return False
    # Constraint 4: Book E cannot be placed in position 2
    if arrangement[1] == 'E':
        return False
    # Constraint 5: Book E cannot be placed in position 3
    if arrangement[2] == 'E':
        return False
    # Constraint 6: Book B must be to the right of book C
    if arrangement.index('B') < arrangement.index('C'):
        return False
    # Constraint 7: Book B cannot be placed in position 3
    if arrangement[2] == 'B':
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break