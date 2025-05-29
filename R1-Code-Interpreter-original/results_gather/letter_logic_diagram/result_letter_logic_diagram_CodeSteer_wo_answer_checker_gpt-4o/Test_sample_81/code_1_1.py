def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Solved the entire grid

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] is not None:
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = None  # Backtrack

    return False

def solve_puzzle():
    grid = [
        [None, 'b', 'd', None, 'e', 'g', 'f'],
        [None, None, 'a', 'e', 'g', None, None],
        ['d', 'a', None, None, None, 'c', None],
        ['a', 'e', None, None, None, 'b', None],
        [None, None, 'f', 'c', 'b', None, 'a'],
        [None, 'f', 'c', None, 'd', None, None],
        ['f', None, 'b', None, None, None, 'g']
    ]

    # Determine the letter for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = {grid[i][j] for i, j in diagonal_indices if grid[i][j] is not None}
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Solve the grid using backtracking
    solve(grid)

    # Format the output
    result = '\n'.join(','.join(str(cell) for cell in row) for row in grid)
    print(f"<<<\n{result}\n>>>")

solve_puzzle()