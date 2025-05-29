def fill_grid(grid):
    # Possible letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_diagonal_letters = letters.copy()
    
    for r, c in minor_diagonal_positions:
        if grid[r][c] != '':
            possible_diagonal_letters.intersection_update(grid[r][c])
    
    # Choose a letter for the minor diagonal
    diagonal_letter = possible_diagonal_letters.pop()
    
    # Fill the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[r])
                col_letters = set(grid[i][c] for i in range(7))
                possible_letters = letters - row_letters - col_letters
                
                # Assign a letter to the cell
                grid[r][c] = possible_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['b', '', '', '', '', 'a', 'f'],
    ['', '', '', '', '', '', 'b'],
    ['', 'e', 'd', 'a', '', '', ''],
    ['e', '', '', 'f', 'b', 'g', ''],
    ['', 'a', '', 'b', 'g', '', 'e'],
    ['', '', 'b', 'g', '', '', 'd'],
    ['', '', '', '', '', 'd', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))