def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['', '', 'f', '', '', 'e', 'c'],
        ['d', 'f', 'b', 'g', 'e', 'c', ''],
        ['f', 'b', 'g', 'e', '', 'a', ''],
        ['', 'g', '', '', 'a', 'd', ''],
        ['', 'e', '', '', 'd', '', 'b'],
        ['e', '', 'a', 'd', 'f', 'b', 'g'],
        ['', '', 'd', 'f', '', 'g', '']
    ]

    # Determine the letter for the minor diagonal
    letters = set('abcdefg')
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_positions if grid[i][j]}
    diagonal_letter = (letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()