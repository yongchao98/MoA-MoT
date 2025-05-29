def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_grid(grid):
    # First, ensure all minor diagonal positions are 'c'
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != 'c':
            return False
        grid[i][j] = 'c'
    
    def solve(row=0, col=0):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True

        # Skip if cell is pre-filled
        if grid[row][col] != '':
            return solve(row, col + 1)

        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    return solve()

# Initialize grid
initial_grid = [
    ['g', 'b', '', '', '', 'a', 'c'],
    ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
    ['', '', '', '', 'c', 'g', ''],
    ['f', '', '', '', '', '', 'e'],
    ['', 'a', 'c', '', 'b', 'e', 'f'],
    ['a', '', 'g', 'b', 'e', 'f', ''],
    ['c', 'g', 'b', '', '', 'd', '']
]

# Create a copy of the grid to work with
grid = [row[:] for row in initial_grid]

# Solve the puzzle
if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution found")