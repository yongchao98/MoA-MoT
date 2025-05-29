import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Check all constraints
    pos = {book: i for i, book in enumerate(arrangement)}
    
    # Constraint 1: A must be to the right of B
    if pos['A'] <= pos['B']:
        return False
    
    # Constraint 2: F must be to the right of A
    if pos['F'] <= pos['A']:
        return False
    
    # Constraint 3: F cannot be in position 3
    if pos['F'] == 2:
        return False
    
    # Constraint 4: A must be adjacent to D
    if abs(pos['A'] - pos['D']) != 1:
        return False
    
    # Constraint 5: C cannot be in position 5
    if pos['C'] == 4:
        return False
    
    # Constraint 6: B must be to the left of F
    if pos['B'] >= pos['F']:
        return False
    
    # Constraint 7: A must be adjacent to C
    if abs(pos['A'] - pos['C']) != 1:
        return False
    
    # Constraint 8: B must be to the left of F (already checked in constraint 6)
    
    # Constraint 9: E cannot be in position 4
    if pos['E'] == 3:
        return False
    
    return True

# Check all permutations
for arrangement in itertools.permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break