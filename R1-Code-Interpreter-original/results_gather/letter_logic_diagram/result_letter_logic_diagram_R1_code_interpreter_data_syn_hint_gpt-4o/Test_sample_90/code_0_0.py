# Initial grid setup
grid = [
    ['c', 'b', '', 'a', 'f', 'g', ''],
    ['b', 'e', 'a', 'f', '', '', ''],
    ['', '', '', 'g', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['', 'g', '', 'c', 'b', '', ''],
    ['g', 'd', '', '', '', '', 'f'],
    ['', '', 'b', '', 'a', '', '']
]

# Function to fill the grid
def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in letters:
        if all(grid[i][n-1-i] == '' or grid[i][n-1-i] == letter for i in range(n)):
            minor_diagonal_letter = letter
            break
    
    # Fill the minor diagonal
    for i in range(n):
        grid[i][n-1-i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                possible_letters = letters - row_letters - col_letters
                grid[i][j] = possible_letters.pop()
    
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))