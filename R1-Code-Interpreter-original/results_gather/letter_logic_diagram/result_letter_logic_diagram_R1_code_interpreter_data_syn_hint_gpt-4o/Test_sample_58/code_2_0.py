def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    minor_diagonal_letter = None
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            minor_diagonal_letter = grid[i][j]
            break
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'  # Default to 'a' if none are pre-filled

    # Try each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if (row, col) in minor_diagonal_indices and letter != minor_diagonal_letter:
            continue  # Skip if the letter doesn't match the minor diagonal requirement

        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

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
solve(grid)

# Format the output
solution = '\n'.join(','.join(row) for row in grid)
print(solution)