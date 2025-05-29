def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)  # Skip pre-filled cells

    # Determine if this cell is on the minor diagonal
    if row + col == 6:
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, minor_diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack

    return False

# Initial grid setup
grid = [
    ['f', '', 'e', '', '', '', ''],
    ['', 'e', '', 'd', '', '', 'f'],
    ['', '', '', '', 'g', 'f', 'a'],
    ['', 'd', 'b', '', '', 'a', ''],
    ['d', 'b', '', '', 'a', 'e', ''],
    ['b', 'g', '', 'a', 'e', 'c', 'd'],
    ['', 'f', 'a', '', '', 'd', '']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    if solve(grid, letter):
        for row in grid:
            print(','.join(row))
        break