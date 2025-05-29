def solve_puzzle(grid):
    from collections import defaultdict

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Update possible values based on the initial grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                letter = grid[r][c]
                possible_values[r][c] = {letter}
                for i in range(7):
                    possible_values[r][i].discard(letter)
                    possible_values[i][c].discard(letter)

    # Function to check if a letter can be placed on the minor diagonal
    def can_place_on_diagonal(letter):
        return all(letter in possible_values[r][6-r] for r in range(7))

    # Function to place a letter and propagate constraints
    def place_letter(r, c, letter):
        grid[r][c] = letter
        possible_values[r][c] = {letter}
        for i in range(7):
            possible_values[r][i].discard(letter)
            possible_values[i][c].discard(letter)

    # Backtracking function with forward checking
    def backtrack():
        # Find the cell with the fewest possibilities
        min_possibilities = 8
        cell = None
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '' and len(possible_values[r][c]) < min_possibilities:
                    min_possibilities = len(possible_values[r][c])
                    cell = (r, c)

        if cell is None:
            return True  # All cells are filled

        r, c = cell
        for letter in possible_values[r][c]:
            # Save the current state
            saved_grid = [row[:] for row in grid]
            saved_possible_values = [row[:] for row in possible_values]

            # Place the letter and propagate
            place_letter(r, c, letter)
            if backtrack():
                return True

            # Restore the state
            grid[:] = saved_grid
            possible_values[:] = saved_possible_values

        return False

    # Try each letter for the minor diagonal
    for diag_letter in 'abcdefg':
        if can_place_on_diagonal(diag_letter):
            # Place diag_letter on the minor diagonal
            for r in range(7):
                place_letter(r, 6-r, diag_letter)

            # Try to fill the rest of the grid
            if backtrack():
                return '\n'.join(','.join(row) for row in grid)

            # Reset the grid and possible values
            grid = [list(row) for row in [
                ['g', '', 'c', 'a', 'd', '', 'e'],
                ['', 'c', '', '', '', '', ''],
                ['', '', '', '', '', 'g', ''],
                ['', '', 'b', '', 'g', 'f', ''],
                ['d', '', '', '', 'f', '', 'a'],
                ['b', '', 'g', 'f', '', 'a', ''],
                ['', '', '', '', 'a', 'd', '']
            ]]
            possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
            for r in range(7):
                for c in range(7):
                    if grid[r][c] != '':
                        letter = grid[r][c]
                        possible_values[r][c] = {letter}
                        for i in range(7):
                            possible_values[r][i].discard(letter)
                            possible_values[i][c].discard(letter)

    return "No solution found"

# Initial grid
grid = [
    ['g', '', 'c', 'a', 'd', '', 'e'],
    ['', 'c', '', '', '', '', ''],
    ['', '', '', '', '', 'g', ''],
    ['', '', 'b', '', 'g', 'f', ''],
    ['d', '', '', '', 'f', '', 'a'],
    ['b', '', 'g', 'f', '', 'a', ''],
    ['', '', '', '', 'a', 'd', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(f"<<<\n{solution}\n>>>")