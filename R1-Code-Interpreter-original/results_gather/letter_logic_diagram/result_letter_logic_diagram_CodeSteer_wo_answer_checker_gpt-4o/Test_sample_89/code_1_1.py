def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid

    # Move to the next cell
    next_row, next_col = (row, col + 1) if col < 6 else (row + 1, 0)

    if grid[row][col] is not None:
        return solve(grid, minor_diagonal_letter, next_row, next_col)

    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if (row, col) in minor_diagonal_positions and letter != minor_diagonal_letter:
            continue
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, minor_diagonal_letter, next_row, next_col):
                return True
            grid[row][col] = None  # Backtrack

    return False

def solve_puzzle():
    grid = [
        [None, None, 'a', 'g', None, 'b', 'c'],
        [None, 'a', None, None, None, 'c', 'd'],
        ['a', 'g', 'e', 'b', None, None, None],
        ['g', None, 'b', None, 'd', 'f', 'a'],
        ['e', 'b', None, None, 'f', None, 'g'],
        [None, None, 'd', 'f', 'a', None, None],
        ['d', 'f', 'a', 'g', 'e', 'b', None]
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    for r, c in minor_diagonal_positions:
        if grid[r][c] in letters:
            letters.remove(grid[r][c])
    minor_diagonal_letter = letters.pop()

    # Fill the minor diagonal with the chosen letter
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter

    # Solve the puzzle using backtracking
    if solve(grid, minor_diagonal_letter):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()