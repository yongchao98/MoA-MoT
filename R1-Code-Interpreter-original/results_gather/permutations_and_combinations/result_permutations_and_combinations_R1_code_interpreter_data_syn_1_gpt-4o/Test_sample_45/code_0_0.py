from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

def is_valid(arrangement):
    # Extract positions of each book
    pos = {book: arrangement.index(book) for book in books}
    
    # Check all constraints
    if pos['B'] <= pos['F']: return False
    if pos['E'] >= pos['F']: return False
    if abs(pos['E'] - pos['C']) != 1: return False
    if abs(pos['B'] - pos['F']) != 1: return False
    if pos['D'] <= pos['E']: return False
    if pos['E'] <= pos['C']: return False
    if pos['C'] >= pos['D']: return False
    
    return True

# Find a valid arrangement
for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break