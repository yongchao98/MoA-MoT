def is_valid(grid, row, col, letter, minor_diagonal_letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    if letter in (grid[i][col] for i in range(7)):
        return False
    # Check minor diagonal constraint
    if (row + col == 6) and (letter != minor_diagonal_letter):
        return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] is None:
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter, minor_diagonal_letter):
                        grid[row][col] = letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = None
                return False
    return True

def solve_puzzle():
    # Initial grid with empty cells represented as None
    grid = [
        ['g', None, None, 'd', None, None, None],
        [None, 'b', 'd', None, None, None, None],
        [None, None, None, None, 'e', 'g', 'a'],
        [None, None, 'f', None, None, 'a', None],
        [None, 'f', 'e', None, None, 'b', None],
        ['f', None, None, None, 'b', 'd', 'c'],
        [None, 'g', 'a', None, None, None, 'f']
    ]

    # Determine the letter for the minor diagonal
    # Choose 'c' as it is not present in any minor diagonal positions
    minor_diagonal_letter = 'c'
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter

    # Solve the puzzle using backtracking
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()