def solve_puzzle(grid):
    import copy
    from itertools import permutations

    # Helper function to check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(row[i] for row in grid)) != 7:  # Check columns
                return False
        return True

    # Helper function to fill the minor diagonal
    def fill_minor_diagonal(grid, letter):
        for i in range(7):
            grid[i][6-i] = letter

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Find the letter for the minor diagonal
    minor_diagonal_indices = [(i, 6-i) for i in range(7)]
    minor_diagonal_letters = {grid[i][6-i] for i in range(7) if grid[i][6-i] != ''}
    all_letters = set('abcdefg')
    minor_diagonal_letter = (all_letters - minor_diagonal_letters).pop()

    # Fill the minor diagonal
    fill_minor_diagonal(grid, minor_diagonal_letter)

    # Fill the rest of the grid
    for i in range(7):
        missing_letters = list(all_letters - set(grid[i]))
        for j in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing_letters.pop()

    # Check if the grid is valid
    if not is_valid(grid):
        return "No solution found"

    # Convert grid back to the required format
    result = '\n'.join(','.join(row) for row in grid)
    return result

# Initial grid
grid = """
,b,,c,d,f,a
,,c,,,,g
,c,,,,g,
c,d,,a,,,e
d,,a,,b,e,
f,a,,,,,
a,,b,e,c,,f
"""

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)