def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter can be placed at grid[row][col]
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    # Check minor diagonal
    if row + col == 6 and letter != diagonal_letter:
        return False
    return True

def solve(grid, row, col, diagonal_letter):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)  # Move to the next row

    if grid[row][col] is not None:
        return solve(grid, row, col + 1, diagonal_letter)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = None  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['b', None, None, None, 'f', None, None],
        [None, 'd', None, None, 'a', None, None],
        ['d', None, None, None, 'g', None, None],
        ['c', 'f', None, 'g', None, 'e', None],
        [None, None, None, None, 'e', 'd', None],
        [None, None, None, 'e', None, 'c', 'f'],
        [None, 'b', 'e', None, 'c', 'f', None]
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] is not None:
            diagonal_letter = grid[i][6-i]
            break

    if diagonal_letter is None:
        diagonal_letter = 'a'  # Default choice if no pre-filled diagonal letter

    # Fill the diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    if solve(grid, 0, 0, diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()