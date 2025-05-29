def fill_grid(grid):
    letters = set('abcdefg')
    n = len(grid)
    
    # Determine the letter for the minor diagonal
    # Check which letter can be placed on the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = letters.copy()
    
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            possible_letters.intersection_update(grid[i][j])
    
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                # Determine the possible letters for this cell
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(n))
                possible_letters = letters - row_letters - col_letters
                
                # Assign the first available letter
                grid[i][j] = possible_letters.pop()
    
    return grid

# Initial grid
grid = [
    ['b', '', '', '', 'f', '', ''],
    ['', 'd', '', '', 'a', '', ''],
    ['d', '', '', '', '', 'g', ''],
    ['c', 'f', '', 'g', '', 'e', ''],
    ['', '', '', '', 'e', 'd', ''],
    ['', '', '', 'e', '', 'c', 'f'],
    ['', 'b', 'e', '', 'c', 'f', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))