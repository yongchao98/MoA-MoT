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

    # Backtracking function to fill the grid
    def backtrack(grid, row, col, minor_diagonal_letter):
        if row == 7:
            return is_valid(grid)

        if col == 7:
            return backtrack(grid, row + 1, 0, minor_diagonal_letter)

        if grid[row][col] != '':
            return backtrack(grid, row, col + 1, minor_diagonal_letter)

        for letter in 'abcdefg':
            if letter in grid[row] or any(grid[i][col] == letter for i in range(7)):
                continue

            if col == 6 - row and letter != minor_diagonal_letter:
                continue

            grid[row][col] = letter
            if backtrack(grid, row, col + 1, minor_diagonal_letter):
                return True
            grid[row][col] = ''

        return False

    # Convert input grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the minor diagonal letter
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

    # Solve the puzzle using backtracking
    if backtrack(grid, 0, 0, minor_diagonal_letter):
        return grid
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