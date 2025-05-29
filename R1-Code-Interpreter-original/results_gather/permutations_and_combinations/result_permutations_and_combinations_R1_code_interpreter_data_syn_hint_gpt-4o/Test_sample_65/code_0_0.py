import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

# Constraints
def is_valid(arrangement):
    pos = {book: arrangement.index(book) for book in books}
    
    # 1. Book A must be to the right of book B.
    if pos['A'] <= pos['B']:
        return False
    
    # 2. Book F must be to the right of book G.
    if pos['F'] <= pos['G']:
        return False
    
    # 3. Book D must be adjacent to book E.
    if abs(pos['D'] - pos['E']) != 1:
        return False
    
    # 4. Book F cannot be placed in position 4.
    if pos['F'] == 3:
        return False
    
    # 5. Book E must be to the left of book A.
    if pos['E'] >= pos['A']:
        return False
    
    # 6. Book C cannot be placed in position 4.
    if pos['C'] == 3:
        return False
    
    # 7. Book B must be to the left of book G.
    if pos['B'] >= pos['G']:
        return False
    
    # 8. Book C must be to the right of book F.
    if pos['C'] <= pos['F']:
        return False
    
    return True

# Generate all permutations of the books
for arrangement in itertools.permutations(books):
    if is_valid(arrangement):
        print(list(arrangement))
        break