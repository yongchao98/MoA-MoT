import itertools

# Books to arrange
books = ['A', 'B', 'C', 'D', 'E', 'G']

# Fixed position for F
fixed_position = 4
fixed_book = 'F'

# Function to check if a given arrangement satisfies all constraints
def is_valid(arrangement):
    # Insert F at the fixed position
    arrangement = arrangement[:fixed_position-1] + [fixed_book] + arrangement[fixed_position-1:]
    
    # Extract positions of each book
    pos = {book: arrangement.index(book) + 1 for book in arrangement}
    
    # Check constraints
    if pos['D'] <= pos['A']:
        return False
    if pos['B'] >= pos['G']:
        return False
    if pos['D'] <= pos['G']:
        return False
    if pos['F'] >= pos['G']:
        return False
    if pos['G'] >= pos['D']:
        return False
    if pos['A'] >= pos['E']:
        return False
    if pos['G'] >= pos['C']:
        return False
    
    return True

# Try all permutations of the remaining books
for perm in itertools.permutations(books):
    if is_valid(list(perm)):
        # Insert F at the fixed position
        result = list(perm[:fixed_position-1]) + [fixed_book] + list(perm[fixed_position-1:])
        print(result)
        break