from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Function to check if a given arrangement satisfies all constraints
def is_valid(arrangement):
    # Constraint 1: Book D must be to the right of book C
    if arrangement.index('D') <= arrangement.index('C'):
        return False
    # Constraint 2: Book G cannot be placed in position 1
    if arrangement[0] == 'G':
        return False
    # Constraint 3: Book F must be to the right of book G
    if arrangement.index('F') <= arrangement.index('G'):
        return False
    # Constraint 4: Book B must be to the right of book G
    if arrangement.index('B') <= arrangement.index('G'):
        return False
    # Constraint 5: Book G must be to the right of book A
    if arrangement.index('G') <= arrangement.index('A'):
        return False
    # Constraint 6: Book D must be placed in position 3
    if arrangement[2] != 'D':
        return False
    # Constraint 7: Book F must be to the right of book E
    if arrangement.index('F') <= arrangement.index('E'):
        return False
    # Constraint 8: Book C cannot be placed in position 1
    if arrangement[0] == 'C':
        return False
    return True

# Find a valid arrangement
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break