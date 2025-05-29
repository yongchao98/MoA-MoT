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

    # Check if the minor diagonal is consistent
    def check_minor_diagonal(grid, letter):
        for i in range(7):
            if grid[i][6-i] != letter:
                return False
        return True

    # Fill the grid
    def fill_grid(grid, row=0):
        if row == 7:
            return is_valid(grid) and check_minor_diagonal(grid, grid[0][6])

        # Find missing letters in the current row
        missing_letters = set('abcdefg') - set(grid[row])
        for perm in permutations(missing_letters):
            new_grid = [list(r) for r in grid]
            idx = 0
            for col in range(7):
                if new_grid[row][col] == '':
                    new_grid[row][col] = perm[idx]
                    idx += 1
            if fill_grid(new_grid, row + 1):
                return new_grid
        return None

    # Start solving from the first row
    solution = fill_grid(grid)
    return solution

# Initial grid
grid = [
    ['','g','','','c','e','a'],
    ['g','','b','','','',''],
    ['','c','e','','',''],
    ['b','','e','a','','g','d'],
    ['','','','g','','b'],
    ['','','','d','b',''],
    ['','','d','b','','']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")