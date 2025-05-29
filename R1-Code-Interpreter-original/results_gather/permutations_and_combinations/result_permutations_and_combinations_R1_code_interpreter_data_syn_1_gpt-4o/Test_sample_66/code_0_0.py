from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    if arrangement.index('G') == 2:  # Constraint 1
        return False
    if arrangement.index('C') > arrangement.index('B'):  # Constraint 2
        return False
    if abs(arrangement.index('B') - arrangement.index('D')) != 1:  # Constraint 3
        return False
    if abs(arrangement.index('G') - arrangement.index('D')) != 1:  # Constraint 4
        return False
    if arrangement.index('D') == 6:  # Constraint 5
        return False
    if arrangement.index('C') > arrangement.index('G'):  # Constraint 6
        return False
    if abs(arrangement.index('G') - arrangement.index('C')) != 1:  # Constraint 7
        return False
    if abs(arrangement.index('B') - arrangement.index('E')) != 1:  # Constraint 8
        return False
    return True

for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break