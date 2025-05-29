def solve_puzzle():
    # Initial grid with empty cells as None
    grid = [
        ['a', 'g', 'b', None, None, None, None],
        [None, None, None, None, 'c', None, None],
        ['b', None, None, None, None, None, 'a'],
        ['f', 'd', 'c', None, None, None, None],
        [None, 'c', 'e', 'a', None, 'b', None],
        [None, None, None, 'g', None, 'f', 'd'],
        ['e', 'a', 'g', None, 'f', None, None]
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = 'c'

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    def is_valid(grid, row, col, letter):
        # Check if the letter is already in the row
        if letter in grid[row]:
            return False
        # Check if the letter is already in the column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] is not None:
            return solve(grid, row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = None

        return False

    solve(grid)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()