def solve_puzzle():
    # Initial grid with pre-filled values
    grid = [
        ['c', 'b', '', '', 'd', '', ''],
        ['b', 'g', 'f', '', 'e', '', 'c'],
        ['g', '', 'd', 'e', 'a', 'c', 'b'],
        ['f', 'd', 'e', '', 'c', 'b', 'g'],
        ['d', 'e', '', '', 'b', 'g', ''],
        ['e', 'a', 'c', 'b', 'g', 'f', 'd'],
        ['', '', 'b', 'g', '', '', 'e']
    ]

    # Determine the letter for the minor diagonal
    letters = set('abcdefg')
    diagonal_letter = None
    for letter in letters:
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Backtracking function to fill the grid
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    solve(grid)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()