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
    if grid[row][col] is not None:  # Skip pre-filled cells
        return solve(grid, row, col + 1)

    # Try each letter from 'a' to 'g'
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            print(f"Trying letter '{letter}' at position ({row}, {col})")
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = None  # Backtrack
            print(f"Backtracking from position ({row}, {col})")

    return False

def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        ['b', None, 'c', 'd', 'f', 'a', None],
        [None, None, 'c', None, None, None, 'g'],
        [None, 'c', None, None, None, 'g', None],
        ['c', 'd', None, 'a', None, None, 'e'],
        ['d', None, 'a', None, 'b', 'e', None],
        ['f', 'a', None, None, None, None, None],
        ['a', None, 'b', 'e', 'c', None, 'f']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    print(f"Using '{diagonal_letter}' for the diagonal")

    # Solve the grid
    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()