def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
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

# Initial grid setup
grid = [
    ['b', '', '', '', 'g', '', 'e'],
    ['', '', 'a', 'g', '', '', ''],
    ['f', 'a', '', '', '', 'b', 'd'],
    ['', 'g', '', '', '', '', ''],
    ['g', '', '', 'b', '', '', ''],
    ['', 'e', 'b', 'd', 'f', 'a', ''],
    ['', 'b', '', '', 'a', '', '']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    if solve(grid, letter):
        break

# Print the solution
for row in grid:
    print(','.join(row))