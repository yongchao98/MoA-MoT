import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Constraint 1: Book B cannot be in position 6
    if arrangement.index('B') == 5:
        return False
    # Constraint 2: Book C cannot be in position 6
    if arrangement.index('C') == 5:
        return False
    # Constraint 3: Book C cannot be in position 5
    if arrangement.index('C') == 4:
        return False
    # Constraint 4: Book B must be adjacent to book A
    if abs(arrangement.index('B') - arrangement.index('A')) != 1:
        return False
    # Constraint 5: Book F must be to the right of book C
    if arrangement.index('F') < arrangement.index('C'):
        return False
    # Constraint 6: Book E must be to the right of book B
    if arrangement.index('E') < arrangement.index('B'):
        return False
    # Constraint 7: Book A must be adjacent to book C
    if abs(arrangement.index('A') - arrangement.index('C')) != 1:
        return False
    # Constraint 8: Book C must be adjacent to book F
    if abs(arrangement.index('C') - arrangement.index('F')) != 1:
        return False
    return True

for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break