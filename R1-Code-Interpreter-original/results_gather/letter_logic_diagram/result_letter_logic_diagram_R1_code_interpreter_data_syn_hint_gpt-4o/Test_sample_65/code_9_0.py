def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the given row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_letter):
    # Try to fill the grid using backtracking
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid, minor_letter):
                            return True
                        grid[row][col] = ''
                return False
    return True

def fill_grid():
    # Initial grid with pre-filled values
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
    for minor_letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(is_valid(grid, i, 6-i, minor_letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = minor_letter
            
            # Try to fill the rest of the grid
            if solve(grid, minor_letter):
                return grid
            
            # Reset the minor diagonal if it doesn't work
            for i in range(7):
                grid[i][6-i] = ''
    return None

# Fill the grid
filled_grid = fill_grid()

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")