def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check diagonal constraint
    if (row + col == 6) and (letter != diagonal_letter):
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
        ['f', None, 'e', None, None, None, None],
        [None, 'e', None, 'd', None, None, 'f'],
        [None, None, None, None, 'g', 'f', 'a'],
        [None, 'd', 'b', None, None, 'a', None],
        ['d', 'b', None, None, 'a', 'e', None],
        ['b', 'g', None, 'a', 'e', 'c', 'd'],
        [None, 'f', 'a', None, None, 'd', None]
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j] is not None)
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    if solve(grid, 0, 0, diagonal_letter):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()