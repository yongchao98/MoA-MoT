from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Extract positions of each book
    pos = {book: arrangement.index(book) for book in books}
    
    # Check all constraints
    if pos['A'] >= pos['B']:  # A must be to the left of B
        return False
    if pos['D'] == 0:  # D cannot be in position 1
        return False
    if pos['C'] == 1:  # C cannot be in position 2
        return False
    if pos['F'] == 1:  # F cannot be in position 2
        return False
    if abs(pos['D'] - pos['A']) != 1:  # D must be adjacent to A
        return False
    if pos['D'] == 5:  # D cannot be in position 6
        return False
    if pos['D'] >= pos['E']:  # D must be to the left of E
        return False
    if pos['C'] == 4:  # C cannot be in position 5
        return False
    
    return True

# Try all permutations and find a valid one
for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break