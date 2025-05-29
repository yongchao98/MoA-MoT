def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
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

def fill_grid():
    grid = [
        ['b', '', '', 'g', 'd', 'c', 'c'],
        ['a', 'f', 'g', '', 'c', 'c', 'e'],
        ['a', 'f', '', '', 'c', 'e', 'b'],
        ['g', 'd', '', 'c', 'e', '', ''],
        ['g', 'd', 'c', 'e', 'b', '', 'f'],
        ['c', 'e', '', 'a', '', '', ''],
        ['', '', 'b', 'a', 'f', '', '']
    ]

    # Pre-fill the minor diagonal with 'c'
    for i in range(7):
        grid[i][6-i] = 'c'

    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

fill_grid()