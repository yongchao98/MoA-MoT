def solve_puzzle(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal_letter = None
    
    # Check if any minor diagonal position is pre-filled
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            minor_diagonal_letter = grid[i][j]
            break
    
    # If no pre-filled letter, choose one that is not in the first row
    if minor_diagonal_letter is None:
        for letter in 'abcdefg':
            if letter not in grid[0]:
                minor_diagonal_letter = letter
                break
    
    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that is not in the current row or column
                for letter in 'abcdefg':
                    if letter not in grid[i] and all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        break
    
    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return result

# Initial grid setup
grid = [
    ['d', '', '', '', '', 'f', ''],
    ['', 'b', '', '', 'f', 'a', 'd'],
    ['b', '', '', 'f', '', '', 'd'],
    ['', '', 'f', '', '', '', ''],
    ['c', '', '', '', 'e', '', 'g'],
    ['', 'a', 'd', '', '', '', 'c'],
    ['', '', 'e', '', '', 'c', 'f']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)