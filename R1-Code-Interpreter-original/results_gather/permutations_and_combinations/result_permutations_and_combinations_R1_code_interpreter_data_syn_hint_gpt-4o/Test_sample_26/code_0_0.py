from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('B') == 5:  # Constraint 1
        return False
    if arrangement.index('C') == 5 or arrangement.index('C') == 4:  # Constraints 2 and 3
        return False
    if abs(arrangement.index('B') - arrangement.index('A')) != 1:  # Constraint 4
        return False
    if arrangement.index('F') < arrangement.index('C'):  # Constraint 5
        return False
    if arrangement.index('E') < arrangement.index('B'):  # Constraint 6
        return False
    if abs(arrangement.index('A') - arrangement.index('C')) != 1:  # Constraint 7
        return False
    if abs(arrangement.index('C') - arrangement.index('F')) != 1:  # Constraint 8
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break