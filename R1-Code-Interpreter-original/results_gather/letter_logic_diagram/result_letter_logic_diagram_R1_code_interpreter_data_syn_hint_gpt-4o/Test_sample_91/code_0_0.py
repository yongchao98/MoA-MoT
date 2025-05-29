def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(n):
        if grid[i][n - 1 - i] != '':
            minor_diagonal_letter = grid[i][n - 1 - i]
            break
    
    # If no pre-filled letter on the minor diagonal, choose one
    if not minor_diagonal_letter:
        minor_diagonal_letter = 'a'  # Arbitrarily choose 'a' if none is pre-filled
    
    # Fill the minor diagonal
    for i in range(n):
        grid[i][n - 1 - i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                possible_letters = letters - row_letters - col_letters
                
                # Assign the first possible letter
                grid[i][j] = possible_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['d', 'g', 'c', 'e', '', 'a', ''],
    ['g', 'c', '', '', '', '', ''],
    ['', '', 'f', '', '', 'd', ''],
    ['e', '', '', '', 'd', 'g', ''],
    ['', '', 'd', 'g', '', '', 'e'],
    ['a', '', '', '', '', '', 'f'],
    ['', '', '', 'e', '', '', 'a']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))