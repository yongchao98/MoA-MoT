def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    if letter in grid[row] or any(grid[r][col] == letter for r in range(7)):
        return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:  # If we've filled all rows, the grid is complete
        return True
    if col == 7:  # Move to the next row
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':  # Skip pre-filled cells
        return solve(grid, row, col + 1)

    # Try placing each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

# Initial grid setup
grid = [
    ['b', '', '', '', 'g', 'e', 'd'],
    ['', '', '', '', '', 'd', 'b'],
    ['', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', 'f'],
    ['', '', 'd', '', '', 'f', ''],
    ['e', 'd', '', '', '', 'f', ''],
    ['d', 'b', '', 'f', '', 'g', 'e']
]

# Choose a letter for the minor diagonal
for letter in 'abcdefg':
    if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
        for r in range(7):
            grid[r][6-r] = letter
        break

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))