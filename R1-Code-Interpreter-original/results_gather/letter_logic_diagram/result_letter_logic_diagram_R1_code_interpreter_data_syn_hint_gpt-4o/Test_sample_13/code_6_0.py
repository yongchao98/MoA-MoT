def is_valid(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid):
    # Choose a letter for the minor diagonal
    for letter in 'abcdefg':
        if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
            for r in range(7):
                grid[r][6-r] = letter
            break

    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    backtrack(0, 0)
    return grid

# Initial grid with pre-filled values
grid = [
    ['a', '', 'd', '', 'e', 'c', ''],
    ['b', 'd', 'f', '', '', '', ''],
    ['d', '', '', 'c', '', '', 'a'],
    ['', '', 'c', 'g', '', '', ''],
    ['', '', 'g', '', 'b', '', ''],
    ['', 'g', '', 'b', '', 'f', ''],
    ['', '', '', 'd', '', 'e', '']
]

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))