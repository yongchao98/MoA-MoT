def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row, col):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    for minor_diagonal_letter in 'abcdefg':
        # Check if the minor diagonal letter is valid
        valid_diagonal = True
        for i in range(7):
            if grid[i][6 - i] != '' and grid[i][6 - i] != minor_diagonal_letter:
                valid_diagonal = False
                break

        if not valid_diagonal:
            continue

        # Set the minor diagonal
        for i in range(7):
            grid[i][6 - i] = minor_diagonal_letter

        # Try to solve the grid with the current minor diagonal letter
        if solve(grid, 0, 0):
            for row in grid:
                print(','.join(row))
            return

    print("No solution found")

# Initial grid with empty cells as ''
grid = [
    ['a', 'e', 'f', 'g', '', 'c', 'd'],
    ['e', 'f', 'g', 'b', '', 'd', 'a'],
    ['', 'g', 'b', '', '', '', 'e'],
    ['', 'b', '', '', 'a', '', 'f'],
    ['', 'c', 'd', '', 'e', '', 'g'],
    ['c', 'd', 'a', 'e', 'f', 'g', 'b'],
    ['d', 'a', 'e', 'f', '', 'b', '']
]

solve_puzzle(grid)