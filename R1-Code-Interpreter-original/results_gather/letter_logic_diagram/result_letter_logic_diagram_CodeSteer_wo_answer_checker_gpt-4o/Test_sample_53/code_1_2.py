def solve_puzzle(grid):
    from copy import deepcopy

    # Function to update possible values based on a new assignment
    def update_possible_values(possible_values, row, col, letter):
        for i in range(7):
            if letter in possible_values[row][i]:
                possible_values[row][i].remove(letter)
            if letter in possible_values[i][col]:
                possible_values[i][col].remove(letter)
        possible_values[row][col] = {letter}

    # Function to enforce the minor diagonal constraint
    def enforce_minor_diagonal(grid, possible_values, letter):
        for i in range(7):
            grid[i][6-i] = letter
            possible_values[i][6-i] = {letter}
            update_possible_values(possible_values, i, 6-i, letter)

    # Backtracking with constraint propagation
    def backtrack(grid, possible_values, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return backtrack(grid, possible_values, row + 1, 0)
        if grid[row][col] != '':
            return backtrack(grid, possible_values, row, col + 1)

        for letter in possible_values[row][col]:
            grid[row][col] = letter
            new_possible_values = deepcopy(possible_values)
            update_possible_values(new_possible_values, row, col, letter)
            if backtrack(grid, new_possible_values, row, col + 1):
                return True
            grid[row][col] = ''
        return False

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Pre-fill the grid and update possible values
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                possible_values[r][c] = {grid[r][c]}
                update_possible_values(possible_values, r, c, grid[r][c])

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        grid_copy = deepcopy(grid)
        possible_values_copy = deepcopy(possible_values)
        enforce_minor_diagonal(grid_copy, possible_values_copy, letter)
        if backtrack(grid_copy, possible_values_copy):
            return grid_copy

    return None

# Initial grid
grid = [
    ['', '', 'f', '', '', '', ''],
    ['', '', 'd', 'e', 'g', '', ''],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', '', 'c', 'b', 'a', ''],
    ['e', '', '', 'b', '', '', 'f'],
    ['', 'c', '', 'a', 'f', '', ''],
    ['', '', '', 'f', 'd', 'e', 'g']
]

solution = solve_puzzle(grid)
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")