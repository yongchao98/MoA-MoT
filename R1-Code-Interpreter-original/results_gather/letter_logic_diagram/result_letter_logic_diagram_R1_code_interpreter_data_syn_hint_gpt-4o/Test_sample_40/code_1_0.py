def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the given row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def find_solution(grid):
    for letter in 'abcdefg':
        # Try each letter for the minor diagonal
        temp_grid = [row[:] for row in grid]
        fill_minor_diagonal(temp_grid, letter)
        if solve(temp_grid):
            return temp_grid
    return None

# Initial grid setup
grid = [
    ['e', '', '', 'g', '', 'c', 'b'],
    ['d', '', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Find the solution
solution = find_solution(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")