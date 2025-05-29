def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
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
    
    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue
        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''
    
    return False

def fill_grid():
    grid = [
        ['', 'e', '', 'g', '', 'c', 'b'],
        ['', 'd', '', '', '', 'b', 'a'],
        ['', '', '', 'c', '', '', ''],
        ['', '', 'c', 'b', '', 'e', 'd'],
        ['f', '', '', '', 'e', 'd', ''],
        ['', '', '', '', '', '', ''],
        ['b', 'a', '', 'd', 'g', '', '']
    ]
    
    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if solve(grid, letter):
            return grid

filled_grid = fill_grid()
for row in filled_grid:
    print(','.join(row))