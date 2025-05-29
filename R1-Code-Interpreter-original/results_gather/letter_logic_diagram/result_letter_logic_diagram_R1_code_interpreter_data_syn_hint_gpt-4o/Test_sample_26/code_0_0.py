def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(i, n - 1 - i) for i in range(n)]
    minor_diagonal_letters = {grid[i][n - 1 - i] for i in range(n) if grid[i][n - 1 - i]}
    minor_diagonal_letter = (letters - minor_diagonal_letters).pop()
    
    # Fill the grid
    for i in range(n):
        row_letters = set(grid[i]) - {''}
        col_letters = set(grid[j][i] for j in range(n)) - {''}
        
        for j in range(n):
            if grid[i][j] == '':
                if (i, j) in minor_diagonal_indices:
                    grid[i][j] = minor_diagonal_letter
                else:
                    possible_letters = letters - row_letters - col_letters
                    grid[i][j] = possible_letters.pop()
                    row_letters.add(grid[i][j])
                    col_letters.add(grid[i][j])
    
    return grid

# Initial grid
grid = [
    ['', 'g', 'a', '', 'd', '', ''],
    ['a', '', '', '', '', 'c', 'f'],
    ['a', '', 'd', 'b', 'c', '', 'g'],
    ['e', 'd', '', 'c', 'f', 'g', 'a'],
    ['', 'b', 'c', '', '', 'a', ''],
    ['', 'c', '', 'g', 'a', 'e', 'd'],
    ['', 'f', '', 'a', '', 'd', 'b']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))