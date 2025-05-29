def solve_puzzle():
    # Initial grid with empty spaces represented by ''
    grid = [
        ['f', '', '', 'c', 'a', '', 'g'],
        ['e', 'b', 'c', 'a', 'd', 'g', ''],
        ['', 'c', '', '', 'g', 'f', 'e'],
        ['c', 'a', '', '', '', 'e', 'b'],
        ['a', 'd', 'g', '', 'e', '', 'c'],
        ['d', 'g', '', 'e', '', 'c', ''],
        ['', '', '', 'b', '', 'a', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][j] != letter for i, j in diagonal_positions):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i, j in diagonal_positions:
        grid[i][j] = diagonal_letter

    # Helper function to check if placing a letter is valid
    def is_valid(row, col, letter):
        # Check row and column constraints
        if letter in grid[row]:
            return False
        if any(grid[i][col] == letter for i in range(7)):
            return False
        return True

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()