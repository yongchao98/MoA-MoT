def solve_puzzle(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal_letter = None
    
    # Check if any pre-filled cell on the minor diagonal can determine the letter
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            minor_diagonal_letter = grid[i][j]
            break
    
    # If no pre-filled letter on the minor diagonal, choose one that doesn't conflict
    if minor_diagonal_letter is None:
        for letter in 'abcdefg':
            if all(grid[i][j] != letter for i, j in minor_diagonal_indices):
                minor_diagonal_letter = letter
                break
    
    # Fill the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # If it's on the minor diagonal, fill with the determined letter
                if (i, j) in minor_diagonal_indices:
                    grid[i][j] = minor_diagonal_letter
                else:
                    # Find a letter that can be placed here
                    used_letters = set(grid[i]) | {grid[k][j] for k in range(7)}
                    for letter in 'abcdefg':
                        if letter not in used_letters:
                            grid[i][j] = letter
                            break
    
    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid with pre-filled values
initial_grid = [
    ['g', '', '', 'd', '', '', ''],
    ['', 'b', 'd', '', '', '', ''],
    ['', '', '', '', 'e', 'g', 'a'],
    ['', '', 'f', '', '', 'a', ''],
    ['', 'f', 'e', '', '', 'b', ''],
    ['f', '', '', '', 'b', 'd', 'c'],
    ['', 'g', 'a', '', '', '', 'f']
]

# Solve the puzzle
solution = solve_puzzle(initial_grid)
print(solution)