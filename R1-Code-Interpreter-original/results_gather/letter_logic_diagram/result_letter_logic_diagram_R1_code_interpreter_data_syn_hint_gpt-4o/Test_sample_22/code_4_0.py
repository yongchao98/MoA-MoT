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
    
    for letter in 'abcdefg':
        if (col == 6 - row and letter == diagonal_letter) or col != 6 - row:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''
    return False

def fill_grid():
    grid = [
        ['c', '', '', '', '', '', ''],
        ['', '', 'e', '', '', '', ''],
        ['d', 'e', 'g', '', 'f', 'c', 'b'],
        ['e', 'g', '', 'f', 'c', '', 'd'],
        ['', 'a', 'f', '', '', '', ''],
        ['', 'f', 'c', 'b', '', 'e', 'g'],
        ['f', 'c', '', '', 'e', 'g', '']
    ]

    # Try each letter for the minor diagonal
    for diagonal_letter in 'abcdefg':
        if all(grid[i][6-i] == '' or grid[i][6-i] == diagonal_letter for i in range(7)):
            if solve(grid, diagonal_letter):
                for row in grid:
                    print(','.join(row))
                return

fill_grid()