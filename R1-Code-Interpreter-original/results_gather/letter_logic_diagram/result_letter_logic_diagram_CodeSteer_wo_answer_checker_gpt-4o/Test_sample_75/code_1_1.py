def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def update_possibilities(grid, possibilities, row, col, letter):
        # Update possibilities for the row and column
        for i in range(7):
            possibilities[row][i].discard(letter)
            possibilities[i][col].discard(letter)

    def solve(grid, possibilities, minor_diagonal_letter):
        # Find the next empty cell with the fewest possibilities
        min_possibilities = 8
        next_cell = None
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '' and (row + col != 6 or grid[row][col] == minor_diagonal_letter):
                    if len(possibilities[row][col]) < min_possibilities:
                        min_possibilities = len(possibilities[row][col])
                        next_cell = (row, col)

        if not next_cell:
            return True  # Solved

        row, col = next_cell
        for letter in possibilities[row][col]:
            if (row + col == 6 and letter != minor_diagonal_letter):
                continue
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                new_possibilities = [row.copy() for row in possibilities]
                update_possibilities(grid, new_possibilities, row, col, letter)
                if solve(grid, new_possibilities, minor_diagonal_letter):
                    return True
                grid[row][col] = ''

        return False

    # Convert input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Initialize possibilities
    possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for row in range(7):
        for col in range(7):
            if grid[row][col] != '':
                possibilities[row][col] = set()

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if solve(grid, possibilities, letter):
            return '\n'.join([','.join(row) for row in grid])

# Input grid
puzzle = """
,,g,d,,,
a,g,,c,,e,
g,,c,f,,,
d,c,f,,,a,
,,,,a,g,d
f,e,,,g,d,
,b,a,g,,c,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(f"<<<\n{solution}\n>>>")