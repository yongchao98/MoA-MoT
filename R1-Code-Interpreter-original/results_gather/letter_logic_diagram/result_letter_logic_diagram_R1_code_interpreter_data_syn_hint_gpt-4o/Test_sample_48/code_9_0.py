def is_valid(grid, row, col, char):
    # Check if char is not in the current row
    if char in grid[row]:
        return False
    # Check if char is not in the current column
    for r in range(7):
        if grid[r][col] == char:
            return False
    return True

def solve(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for char in 'abcdefg':
                    if is_valid(grid, row, col, char):
                        grid[row][col] = char
                        if solve(grid):
                            return True
                        grid[row][col] = ''
                return False
    return True

def fill_grid():
    grid = [
        ['f', 'd', '', '', '', 'c', 'g'],
        ['d', 'e', '', 'a', 'c', '', 'f'],
        ['e', 'b', '', '', 'g', 'f', ''],
        ['', '', 'c', 'g', 'f', '', 'e'],
        ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
        ['', 'g', '', 'd', 'e', '', 'a'],
        ['g', 'f', '', 'e', 'b', 'a', 'c']
    ]

    # Fill the minor diagonal with 'g'
    for i in range(7):
        grid[i][6-i] = 'g'

    if solve(grid):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

fill_grid()