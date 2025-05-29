from itertools import permutations

# Initial grid with pre-filled values
grid = [
    ['','d','','','b','e','g'],
    ['','','','b','','g',''],
    ['f','','','e','g','',''],
    ['a','','','','','d',''],
    ['','g','','','','','a'],
    ['g','','d','','',''],
    ['','c','d','','a','b','']
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to solve the grid
def solve(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0,6), (1,5), (2,4), (3,3), (4,2), (5,1), (6,0)]
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[r][c] in ('', letter) for r, c in minor_diagonal_positions):
            minor_diagonal_letter = letter
            break

    # Fill the minor diagonal
    for r, c in minor_diagonal_positions:
        grid[r][c] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    if backtrack(0, 0):
        return grid
    else:
        return None

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
if solved_grid:
    for row in solved_grid:
        print(','.join(row))
else:
    print("No solution found")