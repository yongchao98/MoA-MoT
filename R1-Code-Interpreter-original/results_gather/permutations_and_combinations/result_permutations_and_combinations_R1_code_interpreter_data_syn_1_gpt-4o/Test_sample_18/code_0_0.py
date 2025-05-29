import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: C cannot be in position 3
    if arrangement[2] == 'C':
        return False
    # Constraint 2: A cannot be in position 2
    if arrangement[1] == 'A':
        return False
    # Constraint 3: C must be to the right of F
    if arrangement.index('C') < arrangement.index('F'):
        return False
    # Constraint 4: F must be to the left of E
    if arrangement.index('F') > arrangement.index('E'):
        return False
    # Constraint 5: A must be to the right of E
    if arrangement.index('A') < arrangement.index('E'):
        return False
    # Constraint 6: F must be to the right of D
    if arrangement.index('F') < arrangement.index('D'):
        return False
    # Constraint 7: A must be adjacent to E
    if abs(arrangement.index('A') - arrangement.index('E')) != 1:
        return False
    return True

for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break