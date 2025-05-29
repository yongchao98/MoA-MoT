def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Try to fill the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                # If it's a minor diagonal cell, it must be the minor_diagonal_letter
                if row + col == 6:
                    grid[row][col] = minor_diagonal_letter
                else:
                    # Try each letter from a to g
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid, minor_diagonal_letter):
                                return True
                            grid[row][col] = ''  # Backtrack
                return False
    return True

# Initial grid setup
grid = [
    ['', 'a', 'e', 'g', 'f', 'b', ''],
    ['a', '', 'g', 'f', '', '', ''],
    ['', '', '', '', 'c', 'd', 'a'],
    ['', 'f', 'b', 'c', '', 'a', ''],
    ['', 'b', 'c', '', 'a', '', 'g'],
    ['b', '', '', '', '', 'g', ''],
    ['c', '', '', 'e', 'g', '', '']
]

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    # Check if the letter is not already in any of the minor diagonal cells
    if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
        # Create a copy of the grid to try this letter
        grid_copy = [row[:] for row in grid]
        if solve(grid_copy, letter):
            # If solved, print the grid
            for row in grid_copy:
                print(','.join(row))
            break