from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: F must be to the right of D
    if arrangement.index('F') <= arrangement.index('D'):
        return False
    # Constraint 2: E must be to the right of D
    if arrangement.index('E') <= arrangement.index('D'):
        return False
    # Constraint 3: B must be adjacent to C
    if abs(arrangement.index('B') - arrangement.index('C')) != 1:
        return False
    # Constraint 4: E cannot be in position 3
    if arrangement.index('E') == 2:
        return False
    # Constraint 5: E must be to the left of F
    if arrangement.index('E') >= arrangement.index('F'):
        return False
    # Constraint 6: E must be adjacent to F
    if abs(arrangement.index('E') - arrangement.index('F')) != 1:
        return False
    # Constraint 7: D must be to the right of A
    if arrangement.index('D') <= arrangement.index('A'):
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break