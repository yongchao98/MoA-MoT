from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    pos = {book: arrangement.index(book) for book in books}
    
    # Check all constraints
    if pos['E'] > pos['F']:
        return False
    if abs(pos['G'] - pos['E']) != 1:
        return False
    if pos['E'] > pos['C']:
        return False
    if pos['D'] > pos['F']:
        return False
    if pos['F'] == 6:
        return False
    if pos['E'] == 6:
        return False
    if pos['A'] < pos['E']:
        return False
    if pos['E'] == 3:
        return False
    if pos['A'] > pos['D']:
        return False
    
    return True

for arrangement in permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break