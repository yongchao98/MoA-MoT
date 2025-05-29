def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve_grid(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve_grid(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve_grid(grid, row, col + 1)

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''

    return False

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        if not is_valid(grid, i, 6 - i, letter):
            return False
        grid[i][6 - i] = letter
    return True

def solve_puzzle(grid):
    for letter in 'abcdefg':
        temp_grid = [row[:] for row in grid]
        if fill_minor_diagonal(temp_grid, letter) and solve_grid(temp_grid):
            return temp_grid
    return None

# Initial grid with empty cells represented as ''
grid = [
    ['f', 'e', 'g', 'a', '', '', ''],
    ['e', 'g', 'a', '', '', '', 'f'],
    ['g', '', '', '', 'b', 'f', 'e'],
    ['', 'd', '', '', 'f', 'e', 'g'],
    ['d', '', 'b', '', 'e', '', 'a'],
    ['', 'b', 'f', 'e', '', 'a', 'd'],
    ['', '', 'e', '', 'a', 'd', 'c']
]

solution = solve_puzzle(grid)

if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")