def is_valid(grid, row, col, char):
    # Check if char is not in the current row
    if char in grid[row]:
        return False
    # Check if char is not in the current column
    for r in range(7):
        if grid[r][col] == char:
            return False
    return True

def solve(grid, minor_diagonal_char):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for char in 'abcdefg':
            if (row + col == 6 and char != minor_diagonal_char) or not is_valid(grid, row, col, char):
                continue
            grid[row][col] = char
            if backtrack(row, col + 1):
                return True
            grid[row][col] = ''
        return False

    return backtrack(0, 0)

# Initial grid setup
grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['a', '', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

# Try each letter for the minor diagonal
for char in 'abcdefg':
    # Copy the grid to try with a new minor diagonal character
    test_grid = [row[:] for row in grid]
    # Set the minor diagonal
    for i in range(7):
        test_grid[i][6-i] = char
    if solve(test_grid, char):
        for row in test_grid:
            print(','.join(row))
        break