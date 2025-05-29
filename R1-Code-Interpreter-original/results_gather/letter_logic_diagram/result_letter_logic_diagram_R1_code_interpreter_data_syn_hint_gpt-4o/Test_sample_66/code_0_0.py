def fill_grid(grid):
    # Define the letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = letters.copy()
    
    for r, c in minor_diagonal_indices:
        if grid[r][c] != '':
            possible_letters.intersection_update(grid[r][c])
    
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for r, c in minor_diagonal_indices:
        grid[r][c] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[r])
                col_letters = set(grid[i][c] for i in range(7))
                possible_letters = letters - row_letters - col_letters
                
                # Assign the first possible letter
                grid[r][c] = possible_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['', 'c', '', 'e', '', '', ''],
    ['c', '', 'e', '', '', '', 'd'],
    ['a', '', 'b', 'f', '', '', 'c'],
    ['', '', '', '', '', 'c', 'a'],
    ['b', '', 'g', 'd', '', '', ''],
    ['', 'g', '', '', '', 'e', 'b'],
    ['', '', 'c', '', 'e', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))