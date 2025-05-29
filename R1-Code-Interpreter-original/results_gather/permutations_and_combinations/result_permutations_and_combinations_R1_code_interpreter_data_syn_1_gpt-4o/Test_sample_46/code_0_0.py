from itertools import permutations

# Fixed positions
fixed_positions = {2: 'F', 4: 'E', 5: 'B'}

# Books to place: A, C, D
books = ['A', 'C', 'D']

# Possible positions for A, C, D
possible_positions = [1, 3, 6]

# Check all permutations of A, C, D in the possible positions
for perm in permutations(possible_positions):
    # Create a shelf with fixed positions
    shelf = [''] * 6
    for pos, book in fixed_positions.items():
        shelf[pos - 1] = book
    
    # Place A, C, D in the current permutation of positions
    for book, pos in zip(books, perm):
        shelf[pos - 1] = book
    
    # Check constraints
    if shelf.index('C') < shelf.index('E') and shelf.index('C') < shelf.index('A') and shelf.index('F') < shelf.index('D'):
        print(shelf)
        break