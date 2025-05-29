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

    # Try placing each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue

        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''  # Backtrack

    return False

# Initial grid setup
grid = [
    ['a', '', 'd', '', 'b', 'g', 'c'],
    ['f', 'd', '', '', '', '', 'a'],
    ['d', '', '', '', '', '', ''],
    ['', '', '', '', '', 'f', ''],
    ['', 'g', '', 'a', 'f', '', 'e'],
    ['', 'c', '', 'f', '', 'e', ''],
    ['', 'a', '', '', 'e', '', '']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    if solve(grid, letter):
        break

# Print the solution
for row in grid:
    print(','.join(row))