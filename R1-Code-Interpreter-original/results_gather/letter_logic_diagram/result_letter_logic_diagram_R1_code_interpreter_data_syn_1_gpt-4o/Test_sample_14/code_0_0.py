def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
                continue
            grid[row][col] = letter
            if backtrack(row, col + 1):
                return True
            grid[row][col] = ''
        return False

    return backtrack(0, 0)

# Initial grid setup
grid = [
    ['', 'f', '', 'c', '', '', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', '', '', '', 'b', ''],
    ['c', 'e', 'g', '', 'b', 'f', ''],
    ['', '', '', '', 'f', 'd', ''],
    ['', 'a', 'b', '', '', 'c', ''],
    ['', 'b', 'f', 'd', '', '', 'g']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    # Copy the grid to avoid modifying the original
    test_grid = [row[:] for row in grid]
    # Set the minor diagonal
    for i in range(7):
        test_grid[i][6-i] = letter
    if solve(test_grid, letter):
        for row in test_grid:
            print(','.join(row))
        break