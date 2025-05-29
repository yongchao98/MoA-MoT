def solve_puzzle(grid):
    from itertools import permutations

    # Check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Fill the grid
    def fill_grid(grid):
        # Find the letter for the minor diagonal
        minor_diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                minor_diagonal_letter = grid[i][6-i]
                break

        # If no letter is found, choose one
        if minor_diagonal_letter is None:
            minor_diagonal_letter = 'a'  # Arbitrarily choose 'a'

        # Fill the minor diagonal
        for i in range(7):
            grid[i][6-i] = minor_diagonal_letter

        # Fill the rest of the grid
        letters = set('abcdefg')
        for i in range(7):
            missing_letters = list(letters - set(grid[i]))
            for j in range(7):
                if grid[i][j] == '':
                    grid[i][j] = missing_letters.pop()

        return grid

    # Convert input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Fill the grid
    filled_grid = fill_grid(grid)

    # Check if the filled grid is valid
    if is_valid(filled_grid):
        return filled_grid
    else:
        return None

# Input grid
input_grid = """
a,g,,,f,,b
,d,,,e,,a
d,,f,e,,a,g
c,f,e,b,a,g,
,,b,,g,d,c
e,b,a,g,d,c,
b,,g,d,c,f,e
"""

# Solve the puzzle
solution = solve_puzzle(input_grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")