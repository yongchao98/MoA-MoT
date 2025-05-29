def solve_puzzle(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal_letter = None
    
    # Check if any minor diagonal position is pre-filled
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            minor_diagonal_letter = grid[i][j]
            break
    
    # If no pre-filled minor diagonal, choose a letter that is not used in any pre-filled cell
    if minor_diagonal_letter is None:
        used_letters = set()
        for row in grid:
            used_letters.update([cell for cell in row if cell != ''])
        for letter in 'abcdefg':
            if letter not in used_letters:
                minor_diagonal_letter = letter
                break
    
    # Fill the minor diagonal with the determined letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(7):
        row_letters = set(grid[i])
        missing_letters = set('abcdefg') - row_letters
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can be placed in this position
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):  # Check column constraint
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break
    
    return grid

# Initial grid setup
grid = [
    ['', 'b', 'd', '', 'e', 'g', 'f'],
    ['', '', 'a', 'e', 'g', '', ''],
    ['d', 'a', '', '', '', 'c', ''],
    ['a', 'e', '', '', '', 'b', ''],
    ['', '', 'f', 'c', 'b', '', 'a'],
    ['', 'f', 'c', '', 'd', '', ''],
    ['f', '', 'b', '', '', '', 'g']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))