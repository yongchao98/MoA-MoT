# Initial grid with empty cells as None
grid = [
    ['d', 'g', 'c', 'e', None, 'a', None],
    ['g', 'c', None, None, None, None, None],
    [None, None, 'f', None, None, 'd', None],
    ['e', None, None, None, 'd', 'g', None],
    [None, None, None, 'd', 'g', None, 'e'],
    ['a', None, None, None, None, None, 'f'],
    [None, None, None, None, 'e', None, 'a']
]

# Letters to be used
letters = {'a', 'b', 'c', 'd', 'e', 'f', 'g'}

# Function to find missing letters in a list
def find_missing(used):
    return list(letters - set(used))

# Fill the grid
for i in range(7):
    # Find missing letters in the row
    row_missing = find_missing(grid[i])
    
    # Find missing letters in the column
    col_missing = find_missing([grid[j][i] for j in range(7)])
    
    # Fill the row and column
    for j in range(7):
        if grid[i][j] is None:
            # Try to fill with a letter that is missing in both row and column
            for letter in row_missing:
                if letter in col_missing:
                    grid[i][j] = letter
                    row_missing.remove(letter)
                    col_missing.remove(letter)
                    break

# Determine the letter for the minor diagonal
# Check which letter can fit all diagonal positions
diagonal_letter = None
for letter in letters:
    if all(grid[i][6-i] is None or grid[i][6-i] == letter for i in range(7)):
        diagonal_letter = letter
        break

# Fill the minor diagonal with the chosen letter
for i in range(7):
    grid[i][6-i] = diagonal_letter

# Print the filled grid
for row in grid:
    print(','.join(row))