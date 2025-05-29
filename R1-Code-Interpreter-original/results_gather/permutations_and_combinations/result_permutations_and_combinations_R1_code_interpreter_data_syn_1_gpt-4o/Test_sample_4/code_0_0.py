from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Extract positions of each book
    pos = {book: arrangement.index(book) for book in books}
    
    # Check all constraints
    if pos['D'] >= pos['F']:  # D must be to the left of F
        return False
    if abs(pos['D'] - pos['A']) != 1:  # D must be adjacent to A
        return False
    if pos['F'] >= pos['C']:  # F must be to the left of C
        return False
    if pos['C'] == 0:  # C cannot be in position 1
        return False
    if pos['B'] == 0:  # B cannot be in position 1
        return False
    if pos['C'] <= pos['B']:  # C must be to the right of B
        return False
    if pos['A'] >= pos['C']:  # A must be to the left of C
        return False
    
    return True

# Check all permutations
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break