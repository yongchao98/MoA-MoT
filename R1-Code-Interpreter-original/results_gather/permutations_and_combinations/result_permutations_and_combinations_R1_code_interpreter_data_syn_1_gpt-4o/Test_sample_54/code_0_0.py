from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Constraint 1: E must be to the right of C
    if arrangement.index('E') <= arrangement.index('C'):
        return False
    # Constraint 2: C cannot be in position 5
    if arrangement.index('C') == 4:
        return False
    # Constraint 3: A cannot be in position 6
    if arrangement.index('A') == 5:
        return False
    # Constraint 4: B cannot be in position 3
    if arrangement.index('B') == 2:
        return False
    # Constraint 5: E must be to the left of D
    if arrangement.index('E') >= arrangement.index('D'):
        return False
    # Constraint 6: D must be to the right of B
    if arrangement.index('D') <= arrangement.index('B'):
        return False
    # Constraint 7: D must be in position 5
    if arrangement.index('D') != 4:
        return False
    # Constraint 8: G must be to the right of F
    if arrangement.index('G') <= arrangement.index('F'):
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break