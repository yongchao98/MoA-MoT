def fill_grid(grid):
    # Define the letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_diagonal_letters = letters.copy()
    
    for r, c in diagonal_positions:
        if grid[r][c] != '':
            possible_diagonal_letters.intersection_update(grid[r][c])
    
    # Choose a letter for the diagonal
    diagonal_letter = possible_diagonal_letters.pop()
    
    # Fill the diagonal with the chosen letter
    for r, c in diagonal_positions:
        grid[r][c] = diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[r])
                col_letters = set(grid[i][c] for i in range(7))
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[r][c] = missing_letter
    
    return grid

# Initial grid
grid = [
    ['', '', 'f', '', '', 'e', 'c'],
    ['d', 'f', 'b', 'g', 'e', 'c', ''],
    ['f', 'b', 'g', 'e', '', 'a', ''],
    ['', 'g', '', '', 'a', 'd', ''],
    ['', 'e', '', '', 'd', '', 'b'],
    ['e', '', 'a', 'd', 'f', 'b', 'g'],
    ['', '', 'd', 'f', '', 'g', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))