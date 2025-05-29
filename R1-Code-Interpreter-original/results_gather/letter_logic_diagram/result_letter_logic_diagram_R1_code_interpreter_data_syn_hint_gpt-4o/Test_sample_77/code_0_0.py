def fill_grid(grid):
    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i][j] in possible_letters:
                possible_letters.remove(grid[i][j])
    
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                col_letters = set('abcdefg') - {grid[k][j] for k in range(7)}
                possible_letters = row_letters & col_letters
                grid[i][j] = possible_letters.pop()
                row_letters.remove(grid[i][j])
    
    return grid

# Initial grid
grid = [
    ['', '', 'd', 'c', 'a', '', ''],
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', '', '', 'g'],
    ['', '', '', '', '', '', 'd'],
    ['', '', '', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))