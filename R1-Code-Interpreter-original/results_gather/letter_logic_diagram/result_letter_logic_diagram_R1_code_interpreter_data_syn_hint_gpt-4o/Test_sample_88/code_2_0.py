from itertools import permutations

# Initial grid setup
grid = [
    ['f', '', '', '', 'e', 'a', 'g'],
    ['c', 'd', 'b', '', 'a', 'g', 'f'],
    ['', 'b', '', 'a', 'g', 'f', 'c'],
    ['b', 'e', '', 'g', 'f', 'c', 'd'],
    ['', 'a', 'g', 'f', 'c', 'd', 'b'],
    ['', 'g', 'f', '', '', 'b', 'e'],
    ['g', '', '', '', '', 'e', 'a']
]

# Function to check if a letter can be placed at grid[row][col]
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to solve the grid using backtracking
def solve_grid(grid):
    # Determine the letter for the minor diagonal
    diagonal_letter = 'f'
    
    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Backtracking function
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
    
    backtrack(0, 0)
    return grid

# Solve the grid
solved_grid = solve_grid(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))