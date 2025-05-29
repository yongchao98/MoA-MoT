from itertools import permutations

# Positions are 1 to 6
positions = [1, 2, 3, 4, 5, 6]

# Fixed positions based on constraints
fixed_positions = {
    'A': 3,
    'D': 4
}

# Books to place: B, C, E, F
books_to_place = ['B', 'C', 'E', 'F']

# Generate all permutations of the books to place
for perm in permutations(books_to_place):
    # Create a list to represent the shelf
    shelf = [None] * 6
    
    # Place fixed books
    shelf[fixed_positions['A'] - 1] = 'A'
    shelf[fixed_positions['D'] - 1] = 'D'
    
    # Place the permutation of books
    j = 0
    for i in range(6):
        if shelf[i] is None:
            shelf[i] = perm[j]
            j += 1
    
    # Check constraints
    if (shelf[1] != 'F' and  # F cannot be in position 2
        shelf.index('C') > shelf.index('D') and  # C must be to the right of D
        abs(shelf.index('F') - shelf.index('B')) == 1 and  # F must be adjacent to B
        shelf.index('B') != 2 and  # B cannot be in position 3
        shelf.index('E') != 3):  # E cannot be in position 4
        print(shelf)
        break