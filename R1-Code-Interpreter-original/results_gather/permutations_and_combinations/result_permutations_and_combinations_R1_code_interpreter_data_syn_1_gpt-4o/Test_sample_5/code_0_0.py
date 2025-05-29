import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']
positions = [1, 2, 3, 4, 5, 6]

# Constraints
def is_valid(arrangement):
    # Constraint 1: Book E cannot be in position 4
    if arrangement.index('E') == 3:
        return False
    # Constraint 2: Book F must be in position 5
    if arrangement.index('F') != 4:
        return False
    # Constraint 3: Book E must be to the left of book D
    if arrangement.index('E') > arrangement.index('D'):
        return False
    # Constraint 4: Book E cannot be in position 6
    if arrangement.index('E') == 5:
        return False
    # Constraint 5: Book A cannot be in position 5
    if arrangement.index('A') == 4:
        return False
    # Constraint 6: Book B cannot be in position 5
    if arrangement.index('B') == 4:
        return False
    # Constraint 7: Book F cannot be in position 2
    if arrangement.index('F') == 1:
        return False
    # Constraint 8: Book A must be to the right of book D
    if arrangement.index('A') < arrangement.index('D'):
        return False
    return True

# Find a valid arrangement
for perm in itertools.permutations(books):
    if is_valid(perm):
        print(list(perm))
        break