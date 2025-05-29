def fill_grid(grid):
    # Possible letters
    letters = set('abcdefg')
    
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters_for_diagonal = letters.copy()
    
    for r, c in minor_diagonal_indices:
        if grid[r][c] != '':
            possible_letters_for_diagonal.intersection_update(grid[r][c])
    
    # Choose a letter for the minor diagonal
    diagonal_letter = possible_letters_for_diagonal.pop()
    
    # Fill the minor diagonal
    for r, c in minor_diagonal_indices:
        grid[r][c] = diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        row_letters = set(grid[r])
        missing_letters = letters - row_letters
        
        for c in range(7):
            if grid[r][c] == '':
                # Find a letter that can be placed here
                column_letters = {grid[i][c] for i in range(7)}
                available_letters = missing_letters - column_letters
                
                # Place the first available letter
                chosen_letter = available_letters.pop()
                grid[r][c] = chosen_letter
                missing_letters.remove(chosen_letter)
    
    return grid

# Initial grid
grid = [
    ['', 'g', '', '', 'c', 'e', 'a'],
    ['g', '', 'b', '', '', '', ''],
    ['', '', 'c', 'e', '', '', ''],
    ['b', '', 'e', 'a', '', 'g', 'd'],
    ['', '', '', '', 'g', '', 'b'],
    ['', '', '', '', 'd', 'b', ''],
    ['', '', '', 'd', 'b', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))