def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    # Try placing each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

# Initial grid with pre-filled letters
grid = [
    ['d', '', '', 'b', 'e', 'g', ''],
    ['', '', 'b', '', 'g', '', ''],
    ['f', '', '', 'e', 'g', '', ''],
    ['a', '', '', '', '', 'd', ''],
    ['', '', 'g', '', '', '', 'a'],
    ['', 'g', '', 'd', '', '', ''],
    ['', 'c', 'd', '', 'a', 'b', '']
]

# Fill the minor diagonal with 'c'
for i in range(7):
    grid[i][6 - i] = 'c'

# Solve the grid
solve(grid)

# Print the completed grid
for row in grid:
    print(','.join(row))