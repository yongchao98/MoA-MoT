def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:  # If we've filled all rows, the grid is complete
        return True
    if col == 7:  # Move to the next row
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':  # Skip pre-filled cells
        return solve(grid, row, col + 1)

    # Try each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def fill_grid():
    # Initial grid with given letters
    grid = [
        ['d', '', '', 'f', 'a', '', ''],
        ['b', '', 'f', '', '', 'e', 'd'],
        ['', 'f', 'a', 'c', 'e', 'd', 'b'],
        ['f', 'a', 'c', '', '', '', ''],
        ['a', '', 'e', '', '', 'g', 'f'],
        ['c', 'e', '', 'b', '', '', 'a'],
        ['e', 'd', '', 'g', 'f', 'a', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = 'c'
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Solve the grid
    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

fill_grid()