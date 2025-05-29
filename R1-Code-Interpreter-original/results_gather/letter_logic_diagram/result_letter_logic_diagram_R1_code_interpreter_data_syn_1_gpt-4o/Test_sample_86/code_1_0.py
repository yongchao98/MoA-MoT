def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in letters:
        if all(grid[i][n-1-i] in ('', letter) for i in range(n)):
            minor_diagonal_letter = letter
            break
    
    # Fill the minor diagonal
    for i in range(n):
        grid[i][n-1-i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                available_letters = letters - row_letters - col_letters
                if available_letters:
                    grid[i][j] = available_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['', 'd', 'c', 'f', '', 'b', 'a'],
    ['', '', 'f', 'e', 'b', 'a', ''],
    ['c', '', 'e', 'b', 'a', '', ''],
    ['', '', 'b', 'a', 'g', 'd', ''],
    ['', '', '', 'g', 'd', 'c', ''],
    ['', '', '', 'd', '', '', 'e'],
    ['', '', '', 'c', 'f', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))