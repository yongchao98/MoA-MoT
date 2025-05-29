def solve_puzzle(grid):
    from copy import deepcopy

    # Initialize possibilities for each cell
    possibilities = [[set('abcdefg') if cell == '' else {cell} for cell in row] for row in grid]

    # Determine the letter for the minor diagonal
    diagonal_letters = [grid[i][6-i] for i in range(7) if grid[i][6-i] != '']
    minor_diagonal_letter = (set('abcdefg') - set(diagonal_letters)).pop()

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
        possibilities[i][6-i] = {minor_diagonal_letter}

    # Update possibilities based on the current grid
    def update_possibilities():
        for i in range(7):
            for j in range(7):
                if len(possibilities[i][j]) == 1:
                    letter = next(iter(possibilities[i][j]))
                    # Remove this letter from the possibilities of other cells in the same row and column
                    for k in range(7):
                        if k != j:
                            possibilities[i][k].discard(letter)
                        if k != i:
                            possibilities[k][j].discard(letter)

    # Check if the current grid satisfies all constraints
    def is_valid():
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Backtracking with constraint propagation
    def backtrack():
        # Find the cell with the fewest possibilities
        min_possibilities = 8
        cell = None
        for i in range(7):
            for j in range(7):
                if len(possibilities[i][j]) > 1 and len(possibilities[i][j]) < min_possibilities:
                    min_possibilities = len(possibilities[i][j])
                    cell = (i, j)

        if cell is None:
            return is_valid()  # All cells are filled

        i, j = cell
        original_possibilities = deepcopy(possibilities[i][j])

        for letter in original_possibilities:
            grid[i][j] = letter
            possibilities[i][j] = {letter}
            update_possibilities()

            if backtrack():
                return True

            # Backtrack
            grid[i][j] = ''
            possibilities[i][j] = original_possibilities
            update_possibilities()

        return False

    update_possibilities()
    if backtrack():
        return grid
    else:
        return None

# Initial grid
grid = [
    ['', '', '', 'b', '', 'e', ''],
    ['f', '', '', 'd', '', 'c', 'a'],
    ['', '', '', 'e', '', 'a', 'f'],
    ['b', 'd', '', '', '', 'f', ''],
    ['d', 'e', '', '', '', '', ''],
    ['', '', '', 'f', 'g', '', ''],
    ['c', '', 'f', '', 'b', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")