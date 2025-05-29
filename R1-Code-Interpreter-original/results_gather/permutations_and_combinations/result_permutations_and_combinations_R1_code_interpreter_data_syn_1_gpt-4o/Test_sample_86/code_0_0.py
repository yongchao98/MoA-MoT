from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

def is_valid(arrangement):
    # Extract positions of each book
    pos = {book: arrangement.index(book) for book in books}
    
    # Check all constraints
    if pos['A'] <= pos['D']:
        return False
    if pos['F'] == 0:
        return False
    if pos['G'] == 4:
        return False
    if abs(pos['B'] - pos['F']) != 1:
        return False
    if pos['B'] == 5:
        return False
    if pos['D'] >= pos['B']:
        return False
    if pos['G'] <= pos['D']:
        return False
    if pos['E'] <= pos['G']:
        return False
    if abs(pos['F'] - pos['A']) != 1:
        return False
    
    return True

# Test all permutations
for perm in permutations(books):
    if is_valid(perm):
        print(list(perm))
        break