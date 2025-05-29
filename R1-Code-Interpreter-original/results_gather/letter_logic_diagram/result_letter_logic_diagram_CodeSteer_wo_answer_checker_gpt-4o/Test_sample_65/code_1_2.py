def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def update_possibilities(grid, possibilities, row, col, letter):
        # Remove the letter from possibilities in the same row and column
        for i in range(7):
            possibilities[row][i].discard(letter)
            possibilities[i][col].discard(letter)

    def solve(grid, possibilities, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, possibilities, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, possibilities, row, col + 1)

        for letter in possibilities[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                new_possibilities = [list(map(set, row)) for row in possibilities]
                update_possibilities(grid, new_possibilities, row, col, letter)
                if solve(grid, new_possibilities, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Convert input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6 - i] != '':
            minor_diagonal_letter = grid[i][6 - i]
            break

    if not minor_diagonal_letter:
        # If no pre-filled letter on the minor diagonal, choose one
        minor_diagonal_letter = 'a'  # Start with 'a' and check if it fits
        for i in range(7):
            if not is_valid(grid, i, 6 - i, minor_diagonal_letter):
                minor_diagonal_letter = chr(ord(minor_diagonal_letter) + 1)

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    # Initialize possibilities
    possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                possibilities[r][c] = set()

    # Solve the rest of the grid
    if not solve(grid, possibilities):
        print("No solution found")

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return f"<<<\n{result}\n>>>"

# Input grid
puzzle = """
a,,d,,b,g,c
f,d,,,,,a
d,,,,,,
,,,,,f,
,g,,a,f,,e
,c,,f,,e,
,a,,,e,,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(solution)