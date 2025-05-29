def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def fill_grid(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row][col] != '':
        return fill_grid(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if can_place(grid, row, col, letter):
            grid[row][col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    return False

def solve_puzzle():
    grid = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            if fill_grid(grid):
                return grid
            for i in range(7):
                grid[i][6-i] = ''

    return None

# Solve the puzzle
solution = solve_puzzle()

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")