def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        # If we've filled all rows, the grid is complete
        return True
    if col == 7:
        # Move to the next row
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] is not None:
        # Skip pre-filled cells
        return solve(grid, minor_diagonal_letter, row, col + 1)

    # Determine the letter to use for the minor diagonal
    if (row + col) == 6:
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = None
    else:
        # Try placing each letter from 'a' to 'g'
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, minor_diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = None

    return False

def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        [None, 'b', None, 'f', 'g', None, 'c'],
        ['b', 'e', None, None, 'a', 'c', None],
        [None, 'f', 'g', 'a', 'c', None, 'b'],
        ['f', None, 'a', None, 'd', None, 'e'],
        ['g', None, None, 'd', None, 'e', 'f'],
        ['a', None, None, None, None, 'f', 'g'],
        ['c', 'd', None, 'e', None, 'g', 'a']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    used_letters = set(grid[i][j] for i, j in minor_diagonal_positions if grid[i][j] is not None)
    all_letters = set('abcdefg')
    minor_diagonal_letter = (all_letters - used_letters).pop()

    # Solve the puzzle
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()