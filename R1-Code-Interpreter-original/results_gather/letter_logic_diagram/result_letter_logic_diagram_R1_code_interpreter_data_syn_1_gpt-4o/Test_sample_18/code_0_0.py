from collections import Counter

# Initial grid setup
grid = [
    ['c', 'b', '', '', 'd', '', ''],
    ['b', 'g', 'f', '', 'e', '', 'c'],
    ['g', '', 'd', 'e', 'a', 'c', 'b'],
    ['f', 'd', 'e', '', 'c', 'b', 'g'],
    ['d', 'e', '', '', 'b', 'g', ''],
    ['e', 'a', 'c', 'b', 'g', 'f', 'd'],
    ['', '', 'b', 'g', '', '', 'e']
]

# Function to find missing letters in a row or column
def find_missing(letters):
    all_letters = set('abcdefg')
    return list(all_letters - set(letters))

# Fill the minor diagonal with a common letter
def fill_minor_diagonal(grid):
    # Find the letter that can fit in all diagonal positions
    possible_letters = set('abcdefg')
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            possible_letters.intersection_update(grid[i][j])
    
    # Choose the first possible letter
    diagonal_letter = possible_letters.pop()
    
    # Fill the diagonal with this letter
    for i in range(7):
        j = 6 - i
        grid[i][j] = diagonal_letter

# Fill the grid
def fill_grid(grid):
    fill_minor_diagonal(grid)
    
    # Fill each row
    for i in range(7):
        missing_letters = find_missing(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing_letters.pop()
    
    # Check columns and adjust if necessary
    for j in range(7):
        column = [grid[i][j] for i in range(7)]
        missing_letters = find_missing(column)
        for i in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing_letters.pop()

fill_grid(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))