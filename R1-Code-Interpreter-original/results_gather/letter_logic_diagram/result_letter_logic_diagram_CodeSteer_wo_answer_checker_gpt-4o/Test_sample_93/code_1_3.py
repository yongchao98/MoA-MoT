def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check the minor diagonal constraint
    if (row + col == 6) and (letter != diagonal_letter):
        return False
    return True

def solve(grid, row=0, col=0, diagonal_letter=None):
    if row == 7:
        # If we've reached the end of the grid, print the solution
        for r in grid:
            print(','.join(r))
        return True

    # Move to the next row if we reach the end of a column
    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)

    # Skip pre-filled cells
    if grid[row][col] is not None:
        return solve(grid, row, col + 1, diagonal_letter)

    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = None  # Backtrack

    return False

def fill_minor_diagonal(grid):
    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    return diagonal_letter

def solve_puzzle():
    grid = [
        ['g', 'b', None, None, None, 'a', 'c'],
        ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
        [None, None, None, None, 'c', 'g', None],
        ['f', None, None, None, None, None, 'e'],
        [None, 'a', 'c', None, 'b', 'e', 'f'],
        ['a', None, 'g', 'b', 'e', 'f', None],
        ['c', 'g', 'b', None, None, 'd', None]
    ]

    # Fill the minor diagonal and get the diagonal letter
    diagonal_letter = fill_minor_diagonal(grid)

    # Solve the puzzle
    if not solve(grid, diagonal_letter=diagonal_letter):
        print("No solution found")

solve_puzzle()