def solve_puzzle(grid):
    from copy import deepcopy

    # Determine possible letters for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_diagonal_letters = set('abcdefg')
    for i, j in diagonal_indices:
        if grid[i][j] != '':
            possible_diagonal_letters.intersection_update(grid[i][j])

    # Function to initialize possible values for each cell
    def initialize_possibilities():
        possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
        for r in range(7):
            for c in range(7):
                if grid[r][c] != '':
                    possibilities[r][c] = set(grid[r][c])
        return possibilities

    # Function to update possibilities after placing a letter
    def update_possibilities(possibilities, row, col, letter):
        for c in range(7):
            possibilities[row][c].discard(letter)
        for r in range(7):
            possibilities[r][col].discard(letter)
        for i, j in diagonal_indices:
            possibilities[i][j].discard(letter)

    # Function to validate the grid
    def validate_grid(grid):
        for r in range(7):
            if set(grid[r]) != set('abcdefg'):
                return False
        for c in range(7):
            if set(grid[r][c] for r in range(7)) != set('abcdefg'):
                return False
        return True

    # Backtracking function with constraint propagation
    def backtrack(grid, possibilities, row, col):
        if row == 7:
            return validate_grid(grid)
        if col == 7:
            return backtrack(grid, possibilities, row + 1, 0)
        if grid[row][col] != '':
            return backtrack(grid, possibilities, row, col + 1)

        for letter in possibilities[row][col]:
            grid[row][col] = letter
            new_possibilities = deepcopy(possibilities)
            update_possibilities(new_possibilities, row, col, letter)
            if backtrack(grid, new_possibilities, row, col + 1):
                return True
            grid[row][col] = ''

        return False

    # Try each possible letter for the diagonal
    for diagonal_letter in possible_diagonal_letters:
        # Fill the diagonal with the current letter
        for i, j in diagonal_indices:
            grid[i][j] = diagonal_letter

        # Initialize possibilities
        possibilities = initialize_possibilities()

        # Start backtracking from the first cell
        if backtrack(grid, possibilities, 0, 0):
            # Format the output
            result = '\n'.join(','.join(row) for row in grid)
            return f"<<<\n{result}\n>>>"

    return "No solution found"

# Initial grid with empty cells as ''
grid = [
    ['', '', 'b', '', '', 'd', 'a'],
    ['', '', 'e', '', 'd', '', 'g'],
    ['', 'e', 'c', 'd', '', '', ''],
    ['', 'c', '', '', '', '', ''],
    ['', 'd', 'a', '', 'f', '', 'e'],
    ['', '', 'g', '', '', '', ''],
    ['', '', 'f', '', 'e', 'c', 'd']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)