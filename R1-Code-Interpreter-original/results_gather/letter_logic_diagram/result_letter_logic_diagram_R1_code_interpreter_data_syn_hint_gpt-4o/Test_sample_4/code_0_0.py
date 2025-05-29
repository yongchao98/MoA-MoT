def fill_grid(grid):
    # Define the letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = letters.copy()
    
    for r, c in minor_diagonal_positions:
        # Remove letters already in the row or column from possible letters
        possible_letters -= set(grid[r]) | {grid[i][c] for i in range(7)}
    
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[r])
                col_letters = {grid[i][c] for i in range(7)}
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[r][c] = missing_letter
    
    return grid

# Initial grid
grid = [
    ['b', '', '', '', 'g', '', 'e'],
    ['', '', 'a', 'g', '', '', ''],
    ['f', 'a', '', '', '', 'b', 'd'],
    ['', 'g', '', '', '', '', ''],
    ['g', '', '', 'b', '', '', ''],
    ['', 'e', 'b', 'd', 'f', 'a', ''],
    ['', 'b', '', '', 'a', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))