def fill_grid(grid):
    # Define the letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = letters.copy()
    
    for r, c in minor_diagonal_positions:
        row_letters = set(grid[r])
        col_letters = set(grid[i][c] for i in range(7))
        possible_letters &= (letters - row_letters - col_letters)
    
    # Choose the letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                row_letters = set(grid[r])
                col_letters = set(grid[i][c] for i in range(7))
                available_letters = letters - row_letters - col_letters
                grid[r][c] = available_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))