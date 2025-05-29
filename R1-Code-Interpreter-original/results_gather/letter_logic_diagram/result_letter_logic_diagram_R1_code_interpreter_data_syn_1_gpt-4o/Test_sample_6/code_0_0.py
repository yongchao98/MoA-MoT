def solve_puzzle(grid):
    import copy
    
    # Helper function to check if a letter can be placed at a given position
    def can_place(grid, row, col, letter):
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
            if all(can_place(grid, r, 6-r, letter) for r in range(7)):
                return letter
        return None

    # Fill the grid
    def fill_grid(grid, minor_diagonal_letter):
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    for letter in 'abcdefg':
                        if can_place(grid, r, c, letter):
                            grid[r][c] = letter
                            break
        return grid

    # Copy the grid to avoid modifying the original
    grid = copy.deepcopy(grid)
    
    # Find the letter for the minor diagonal
    minor_diagonal_letter = find_minor_diagonal_letter(grid)
    
    # Fill the minor diagonal
    for r in range(7):
        grid[r][6-r] = minor_diagonal_letter
    
    # Fill the rest of the grid
    filled_grid = fill_grid(grid, minor_diagonal_letter)
    
    return filled_grid

# Initial grid
grid = [
    ['c', 'b', 'a', 'd', '', 'f', ''],
    ['b', '', '', 'e', '', '', ''],
    ['', 'd', '', 'f', 'g', '', ''],
    ['d', '', '', 'g', '', '', ''],
    ['', '', 'g', '', '', '', ''],
    ['', 'g', 'c', '', 'a', '', ''],
    ['', '', 'b', '', 'd', '', 'f']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
for row in solution:
    print(','.join(row))