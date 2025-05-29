from collections import Counter

# Initial grid setup
grid = [
    ['a', 'e', 'f', 'g', '', 'c', 'd'],
    ['e', 'f', 'g', 'b', '', 'd', 'a'],
    ['', 'g', 'b', '', '', '', 'e'],
    ['', 'b', '', '', 'a', '', 'f'],
    ['', 'c', 'd', '', 'e', '', 'g'],
    ['c', 'd', 'a', 'e', 'f', 'g', 'b'],
    ['d', 'a', 'e', 'f', '', 'b', '']
]

# Function to find missing letters in a list
def find_missing_letters(lst):
    all_letters = set('abcdefg')
    return list(all_letters - set(lst))

# Fill the minor diagonal with a consistent letter
# Check which letter can be used for the minor diagonal
minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
possible_letters = set('abcdefg')

# Remove letters already in the minor diagonal positions
for i, j in minor_diagonal_indices:
    if grid[i][j]:
        possible_letters.discard(grid[i][j])

# Choose a letter for the minor diagonal
minor_diagonal_letter = possible_letters.pop()

# Fill the minor diagonal
for i, j in minor_diagonal_indices:
    grid[i][j] = minor_diagonal_letter

# Fill the rest of the grid
for i in range(7):
    for j in range(7):
        if grid[i][j] == '':
            # Find missing letters in the row and column
            row_letters = set(grid[i])
            col_letters = set(grid[k][j] for k in range(7))
            missing_letters = find_missing_letters(row_letters | col_letters)
            grid[i][j] = missing_letters[0]  # Choose the first available letter

# Print the completed grid
for row in grid:
    print(','.join(row))