def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Solved

    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue

        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''

    return False

def find_minor_diagonal_letter(grid):
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            return grid[i][j]
    return 'a'  # Default to 'a' if no pre-filled letter

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

# Determine the minor diagonal letter
minor_diagonal_letter = find_minor_diagonal_letter(grid)

# Solve the puzzle
solve(grid, minor_diagonal_letter)

# Print the solved grid
for row in grid:
    print(','.join(row))