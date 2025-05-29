def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def update_possibilities(grid, possibilities, row, col, letter, add=True):
    # Update the possibilities for the row and column
    for i in range(7):
        if add:
            possibilities[row][i].discard(letter)
            possibilities[i][col].discard(letter)
        else:
            possibilities[row][i].add(letter)
            possibilities[i][col].add(letter)

def solve(grid, possibilities, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        # If we've filled all rows, the grid is complete
        return True
    if col == 7:
        # Move to the next row
        return solve(grid, possibilities, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] is not None:
        # Skip pre-filled cells
        return solve(grid, possibilities, minor_diagonal_letter, row, col + 1)

    # Determine the letter to use for the minor diagonal
    if (row + col) == 6:
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            update_possibilities(grid, possibilities, row, col, minor_diagonal_letter)
            if solve(grid, possibilities, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = None
            update_possibilities(grid, possibilities, row, col, minor_diagonal_letter, add=False)
    else:
        # Try placing each possible letter
        for letter in possibilities[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                update_possibilities(grid, possibilities, row, col, letter)
                if solve(grid, possibilities, minor_diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = None
                update_possibilities(grid, possibilities, row, col, letter, add=False)

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

    # Initialize possibilities for each cell
    possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if grid[i][j] is not None:
                update_possibilities(grid, possibilities, i, j, grid[i][j])

    # Solve the puzzle
    if solve(grid, possibilities, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()