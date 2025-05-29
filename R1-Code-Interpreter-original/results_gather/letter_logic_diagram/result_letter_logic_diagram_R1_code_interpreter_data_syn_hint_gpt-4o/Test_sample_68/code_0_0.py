def fill_grid(grid):
    # Possible letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for letter in letters:
        if all(grid[r][c] == '' or grid[r][c] == letter for r, c in diagonal_positions):
            diagonal_letter = letter
            break
    
    # Fill the minor diagonal with the chosen letter
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
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))