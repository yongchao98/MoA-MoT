def solve_puzzle(grid):
    # Determine the letters that need to be on the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    all_letters = set('abcdefg')
    
    # Find the letter that can be used for the minor diagonal
    possible_diagonal_letters = all_letters.copy()
    for r, c in minor_diagonal_indices:
        if grid[r][c] != '':
            possible_diagonal_letters.intersection_update(grid[r][c])
    
    # Choose one letter for the minor diagonal
    diagonal_letter = possible_diagonal_letters.pop()
    
    # Fill the minor diagonal with the chosen letter
    for r, c in minor_diagonal_indices:
        grid[r][c] = diagonal_letter
    
    # Fill the rest of the grid
    for r in range(7):
        row_letters = set(grid[r])
        missing_row_letters = all_letters - row_letters
        
        for c in range(7):
            if grid[r][c] == '':
                # Find the missing letter for this column
                col_letters = set(grid[i][c] for i in range(7))
                missing_col_letters = all_letters - col_letters
                
                # Find the intersection of missing row and column letters
                possible_letters = missing_row_letters.intersection(missing_col_letters)
                
                # Fill the cell with one of the possible letters
                grid[r][c] = possible_letters.pop()
    
    return grid

# Initial grid with empty cells as ''
grid = [
    ['d', 'a', '', '', 'g', '', ''],
    ['', '', '', 'g', 'e', '', ''],
    ['f', '', 'g', '', 'c', '', ''],
    ['', 'g', '', '', '', 'a', ''],
    ['g', 'e', '', '', '', 'f', ''],
    ['e', '', '', '', '', 'b', ''],
    ['', 'd', '', 'f', 'b', '', 'e']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))