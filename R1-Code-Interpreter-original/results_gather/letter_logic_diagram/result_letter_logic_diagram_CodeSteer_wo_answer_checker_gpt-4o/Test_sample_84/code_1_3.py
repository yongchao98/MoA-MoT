def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, possible_values, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, possible_values, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, possible_values, row, col + 1)  # Skip pre-filled cells

    for letter in possible_values[row][col]:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            # Update possible values for the row and column
            original_row_values = [set(possible_values[row][i]) for i in range(7)]
            original_col_values = [set(possible_values[i][col]) for i in range(7)]
            for i in range(7):
                possible_values[row][i].discard(letter)
                possible_values[i][col].discard(letter)
            if solve(grid, possible_values, row, col + 1):
                return True
            # Backtrack
            grid[row][col] = ''
            for i in range(7):
                possible_values[row][i] = original_row_values[i]
                possible_values[i][col] = original_col_values[i]

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

    # Initialize possible values for each cell
    possible_values = [[set(letters) for _ in range(7)] for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                possible_values[i][j] = set()

    if solve(grid, possible_values):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()