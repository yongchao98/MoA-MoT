def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        # If we've filled all rows, the grid is complete
        return True
    if col == 7:
        # Move to the next row
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        # Skip pre-filled cells
        return solve(grid, row, col + 1)

    # Determine the letter for the minor diagonal
    if row + col == 6:
        # If on the minor diagonal, use the diagonal letter
        diagonal_letter = grid[0][6]
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
        return False

    # Try placing each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''

    return False

def solve_puzzle():
    grid = [
        ['c', '', '', '', '', '', ''],
        ['', '', 'e', '', '', '', ''],
        ['d', 'e', 'g', '', 'f', 'c', 'b'],
        ['e', 'g', '', 'f', 'c', '', 'd'],
        ['', 'a', 'f', '', '', '', ''],
        ['', 'f', 'c', 'b', '', 'e', 'g'],
        ['f', 'c', '', '', 'e', 'g', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the first cell of the minor diagonal with the chosen letter
    grid[0][6] = diagonal_letter

    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()