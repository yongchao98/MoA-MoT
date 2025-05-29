def solve_puzzle(grid):
    from copy import deepcopy

    # Helper function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Find the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    used_letters = set(grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != '')
    all_letters = set('abcdefg')
    minor_diagonal_letter = (all_letters - used_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter

    # Initialize available letters for each row and column
    row_available = [set('abcdefg') - set(row) for row in grid]
    col_available = [set('abcdefg') - set(grid[i][j] for i in range(7)) for j in range(7)]

    # Try to fill the rest of the grid
    def backtrack(row=0, col=0):
        if row == 7:
            return is_valid(grid)
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in row_available[row] & col_available[col]:
            grid[row][col] = letter
            row_available[row].remove(letter)
            col_available[col].remove(letter)

            if backtrack(row, col + 1):
                return True

            # Backtrack
            grid[row][col] = ''
            row_available[row].add(letter)
            col_available[col].add(letter)

        return False

    # Start backtracking from the first cell
    backtrack()

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid setup
grid = [
    ['', 'g', '', '', 'c', 'e', 'a'],
    ['g', '', 'b', '', '', '', ''],
    ['', '', 'c', 'e', '', '', ''],
    ['b', '', 'e', 'a', '', 'g', 'd'],
    ['', '', '', '', 'g', '', 'b'],
    ['', '', '', '', 'd', 'b', ''],
    ['', '', '', 'd', 'b', '', '']
]

# Solve the puzzle
solution = solve_puzzle(deepcopy(grid))
print("<<<")
print(solution)
print(">>>")