from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: B must be adjacent to C
    if not (arrangement.index('B') == arrangement.index('C') - 1 or arrangement.index('C') == arrangement.index('B') - 1):
        return False
    # Constraint 2: A must be to the left of C
    if arrangement.index('A') > arrangement.index('C'):
        return False
    # Constraint 3: E cannot be in position 3
    if arrangement[2] == 'E':
        return False
    # Constraint 4: F must be to the right of B
    if arrangement.index('F') < arrangement.index('B'):
        return False
    # Constraint 5: F must be to the right of C
    if arrangement.index('F') < arrangement.index('C'):
        return False
    # Constraint 6: F cannot be in position 2
    if arrangement[1] == 'F':
        return False
    # Constraint 7: A cannot be in position 4
    if arrangement[3] == 'A':
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break