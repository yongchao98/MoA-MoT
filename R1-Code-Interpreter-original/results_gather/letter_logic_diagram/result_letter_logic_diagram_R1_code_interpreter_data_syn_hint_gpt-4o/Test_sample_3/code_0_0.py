def solve_puzzle(grid):
    import copy
    
    # Helper function to check if placing a letter is valid
    def is_valid(grid, row, col, letter):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    # Try to find a suitable letter for the minor diagonal
    def find_minor_diagonal_letter(grid):
        for letter in 'abcdefg':
            valid = True
            for i in range(7):
                if not is_valid(grid, i, 6-i, letter):
                    valid = False
                    break
            if valid:
                return letter
        return None

    # Fill the grid
    def fill_grid(grid, minor_diagonal_letter):
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    for letter in 'abcdefg':
                        if is_valid(grid, i, j, letter):
                            grid[i][j] = letter
                            break
        return grid

    # Copy the grid to avoid modifying the original
    grid_copy = copy.deepcopy(grid)
    
    # Find the letter for the minor diagonal
    minor_diagonal_letter = find_minor_diagonal_letter(grid_copy)
    
    # Fill the minor diagonal
    for i in range(7):
        grid_copy[i][6-i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    filled_grid = fill_grid(grid_copy, minor_diagonal_letter)
    
    return filled_grid

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

# Print the solution
for row in solution:
    print(','.join(row))