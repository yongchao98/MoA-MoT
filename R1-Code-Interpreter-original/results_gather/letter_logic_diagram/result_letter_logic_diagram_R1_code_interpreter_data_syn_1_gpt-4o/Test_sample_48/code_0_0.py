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

    # Check if the minor diagonal has the same letter
    def check_minor_diagonal(grid, letter):
        return all(grid[i][6-i] == letter for i in range(7))

    # Try to fill the grid
    def fill_grid(grid, row=0):
        if row == 7:
            # Check if the minor diagonal is valid
            if check_minor_diagonal(grid, grid[0][6]):
                return grid
            return None

        # Find missing letters in the current row
        missing_letters = set('abcdefg') - set(grid[row])
        for perm in permutations(missing_letters):
            new_row = list(grid[row])
            perm_index = 0
            for col in range(7):
                if new_row[col] == '_':
                    new_row[col] = perm[perm_index]
                    perm_index += 1
            grid[row] = new_row
            if is_valid(grid):
                result = fill_grid(grid, row + 1)
                if result:
                    return result
            grid[row] = list(grid[row])  # Reset row
        return None

    # Convert input grid to a list of lists
    grid = [list(row) for row in grid]
    solution = fill_grid(grid)
    return solution

# Initial grid with underscores for empty cells
initial_grid = [
    ['f', 'd', '_', '_', '_', 'c', 'g'],
    ['d', 'e', '_', 'a', 'c', '_', 'f'],
    ['e', 'b', '_', '_', 'g', 'f', '_'],
    ['_', '_', 'c', 'g', 'f', '_', 'e'],
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['_', 'g', '_', 'd', 'e', '_', 'a'],
    ['g', 'f', '_', 'e', 'b', 'a', 'c']
]

solution = solve_puzzle(initial_grid)
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")