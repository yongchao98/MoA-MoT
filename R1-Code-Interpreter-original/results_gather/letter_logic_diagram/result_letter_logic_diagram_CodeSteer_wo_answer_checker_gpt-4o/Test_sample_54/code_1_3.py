def is_valid(grid, row, col, letter, minor_diagonal_letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check the minor diagonal condition
    if col == 6 - row and letter != minor_diagonal_letter:
        return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        # If we've reached the end of the grid, print the solution
        for r in grid:
            print(','.join(r))
        return True

    # Move to the next row if we reach the end of a column
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)

    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, minor_diagonal_letter):
            # Debug statement to trace the execution
            print(f"Trying letter '{letter}' at position ({row}, {col})")

            # Place the letter and move to the next cell
            grid[row][col] = letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True

            # Backtrack
            grid[row][col] = ''
            print(f"Backtracking from position ({row}, {col})")

    return False

def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['', '', '', '', '', 'd', ''],
        ['g', '', '', '', 'd', 'f', ''],
        ['', '', '', '', 'f', 'e', 'a'],
        ['b', '', '', 'f', 'e', '', ''],
        ['c', '', '', '', 'a', 'g', ''],
        ['d', 'f', '', '', 'g', '', ''],
        ['', 'e', 'a', '', 'b', '', '']
    ]

    # Determine the letter for the minor diagonal
    # Based on the pre-filled grid, 'f' is consistent with the minor diagonal
    minor_diagonal_letter = 'f'

    # Solve the puzzle
    if not solve(grid, minor_diagonal_letter):
        print("No solution found.")

solve_puzzle()