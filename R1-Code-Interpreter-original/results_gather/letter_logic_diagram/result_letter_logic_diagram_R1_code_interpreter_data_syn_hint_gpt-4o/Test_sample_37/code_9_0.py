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
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                # Try placing each character from 'a' to 'g'
                for char in 'abcdefg':
                    if is_valid(grid, row, col, char):
                        # Place the character
                        grid[row][col] = char
                        # Check if this is a minor diagonal cell
                        if row + col == 6:
                            if char != minor_diagonal_char:
                                grid[row][col] = ''
                                continue
                        # Recursively try to fill the rest of the grid
                        if solve(grid, minor_diagonal_char):
                            return True
                        # Backtrack
                        grid[row][col] = ''
                return False
    return True

# Initial grid setup
grid = [
    ['c', '', 'g', '', 'f', 'e', ''],
    ['a', 'g', '', '', '', 'd', ''],
    ['', '', '', '', 'd', '', ''],
    ['', '', '', 'd', '', 'a', 'g'],
    ['', 'e', '', 'c', 'a', 'g', 'b'],
    ['e', '', '', 'a', 'g', 'b', 'f'],
    ['d', '', 'a', '', 'b', '', '']
]

# Try each letter for the minor diagonal
for char in 'abcdefg':
    # Copy the grid to try with a new minor diagonal character
    test_grid = [row[:] for row in grid]
    # Set the minor diagonal with the current character
    for i in range(7):
        test_grid[i][6-i] = char
    if solve(test_grid, char):
        # If solved, print the grid
        for row in test_grid:
            print(','.join(row))
        break