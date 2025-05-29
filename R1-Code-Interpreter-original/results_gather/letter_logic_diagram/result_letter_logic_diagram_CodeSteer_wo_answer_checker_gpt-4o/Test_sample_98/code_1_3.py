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

    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def fill_grid():
    # Initial grid with empty cells represented as ''
    grid = [
        ['g', 'c', 'b', 'f', 'e', 'd', ''],
        ['', '', '', 'e', 'd', 'a', ''],
        ['', '', 'e', 'd', 'a', 'g', ''],
        ['f', 'e', '', '', '', 'c', 'b'],
        ['e', 'd', 'a', '', '', '', 'f'],
        ['d', '', 'g', 'c', '', '', 'e'],
        ['', 'g', '', 'b', 'f', '', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j])
    all_letters = set('abcdefg')
    missing_diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the missing letter
    for i, j in diagonal_positions:
        grid[i][j] = missing_diagonal_letter

    # Solve the grid
    if solve(grid):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

fill_grid()