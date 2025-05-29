def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, diagonal_letter, row, col + 1)
    
    if row + col == 6:  # Minor diagonal condition
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''
    return False

def fill_grid():
    grid = [
        ['', '', 'a', 'g', '', 'b', 'c'],
        ['', 'a', '', '', '', 'c', 'd'],
        ['a', 'g', 'e', 'b', '', '', ''],
        ['g', '', 'b', '', 'd', 'f', 'a'],
        ['e', 'b', '', '', 'f', '', 'g'],
        ['', '', 'd', 'f', 'a', '', ''],
        ['', 'd', 'f', 'a', 'g', 'e', 'b']
    ]

    for diagonal_letter in 'abcdefg':
        if solve(grid, diagonal_letter):
            for row in grid:
                print(','.join(row))
            return

fill_grid()