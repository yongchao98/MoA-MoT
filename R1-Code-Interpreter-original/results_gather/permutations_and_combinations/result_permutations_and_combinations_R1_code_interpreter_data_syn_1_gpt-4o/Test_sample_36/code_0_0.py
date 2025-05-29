from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: E must be to the right of A
    if arrangement.index('E') < arrangement.index('A'):
        return False
    # Constraint 2: A cannot be in position 6
    if arrangement.index('A') == 5:
        return False
    # Constraint 3: C must be adjacent to D
    if abs(arrangement.index('C') - arrangement.index('D')) != 1:
        return False
    # Constraint 4: E must be to the left of C
    if arrangement.index('E') > arrangement.index('C'):
        return False
    # Constraint 5: D cannot be in position 6
    if arrangement.index('D') == 5:
        return False
    # Constraint 6: E must be adjacent to D
    if abs(arrangement.index('E') - arrangement.index('D')) != 1:
        return False
    # Constraint 7: E must be adjacent to F
    if abs(arrangement.index('E') - arrangement.index('F')) != 1:
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break