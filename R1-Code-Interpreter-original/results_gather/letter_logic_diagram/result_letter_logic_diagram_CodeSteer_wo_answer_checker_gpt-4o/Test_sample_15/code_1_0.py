def solve_puzzle(grid):
    from collections import defaultdict

    def initialize_possibilities(grid):
        possibilities = defaultdict(lambda: set('abcdefg'))
        for r in range(7):
            for c in range(7):
                if grid[r][c] != '':
                    possibilities[(r, c)] = {grid[r][c]}
                    for i in range(7):
                        possibilities[(r, i)].discard(grid[r][c])
                        possibilities[(i, c)].discard(grid[r][c])
        return possibilities

    def is_valid(grid, row, col, letter):
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def solve(grid, possibilities, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, possibilities, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, possibilities, row, col + 1)

        for letter in possibilities[(row, col)]:
            if (row + col == 6 and letter != diagonal_letter) or not is_valid(grid, row, col, letter):
                continue
            grid[row][col] = letter
            new_possibilities = {k: v.copy() for k, v in possibilities.items()}
            for i in range(7):
                new_possibilities[(row, i)].discard(letter)
                new_possibilities[(i, col)].discard(letter)
            if solve(grid, new_possibilities, row, col + 1):
                return True
            grid[row][col] = ''
        return False

    # Convert input grid to a list of lists
    grid = [list(row.split(',')) for row in grid.strip().split('\n')]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6 - i] in ('', letter) for i in range(7)):
            diagonal_letter = letter
            break

    # Initialize possibilities
    possibilities = initialize_possibilities(grid)

    # Solve the puzzle
    if solve(grid, possibilities):
        # Format the output
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

# Input grid
puzzle = """
c,,,,a,,b
g,e,f,a,,,
,,,d,,,g
f,a,,,c,g,
,d,,c,,e,f
d,,,,,f,
,,,,,,
"""

solve_puzzle(puzzle)