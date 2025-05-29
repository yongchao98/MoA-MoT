def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the current row and column
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
    
    if (row, col) in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
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

def solve_puzzle():
    grid = [
        ['d', 'f', 'e', 'b', 'a', 'g', ''],
        ['f', '', 'b', 'a', 'g', '', ''],
        ['e', 'b', '', '', 'c', 'd', ''],
        ['b', '', 'g', 'c', '', 'f', ''],
        ['b', 'a', '', 'c', 'd', '', ''],
        ['a', '', '', 'd', 'f', 'e', 'b'],
        ['', '', '', '', 'e', '', 'a']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    for r, c in minor_diagonal_positions:
        if grid[r][c] != '':
            letters.discard(grid[r][c])
    minor_diagonal_letter = letters.pop()

    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()