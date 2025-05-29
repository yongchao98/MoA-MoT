def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)
    
    if col == 6 - row:  # Minor diagonal
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, minor_diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''
    return False

def fill_grid():
    grid = [
        ['d', '', '', 'e', 'b', '', 'g'],
        ['', '', 'e', 'b', '', '', 'd'],
        ['', '', 'b', '', '', '', 'c'],
        ['', '', 'a', '', 'd', 'c', 'f'],
        ['', 'a', 'g', 'd', '', '', 'e'],
        ['', '', '', '', '', '', ''],
        ['', '', '', 'f', '', '', 'a']
    ]

    for minor_diagonal_letter in 'abcdefg':
        if solve(grid, minor_diagonal_letter):
            for row in grid:
                print(','.join(row))
            return

fill_grid()