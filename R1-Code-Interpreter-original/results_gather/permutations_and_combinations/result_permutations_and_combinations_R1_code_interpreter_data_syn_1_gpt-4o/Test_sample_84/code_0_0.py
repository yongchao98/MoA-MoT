from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Extract positions of each book
    pos = {book: arrangement.index(book) for book in books}
    
    # Check all constraints
    if pos['B'] >= pos['G']:
        return False
    if pos['B'] >= pos['E']:
        return False
    if abs(pos['F'] - pos['B']) != 1:
        return False
    if abs(pos['G'] - pos['C']) != 1:
        return False
    if pos['C'] >= pos['D']:
        return False
    if pos['C'] == 4:  # position 5 is index 4
        return False
    if pos['C'] >= pos['A']:
        return False
    if pos['A'] == 1:  # position 2 is index 1
        return False
    
    return True

# Try all permutations and find a valid one
for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break