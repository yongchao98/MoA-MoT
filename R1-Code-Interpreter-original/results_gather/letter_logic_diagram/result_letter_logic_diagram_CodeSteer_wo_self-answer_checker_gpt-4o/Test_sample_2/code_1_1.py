def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter, diagonal_positions, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, diagonal_letter, diagonal_positions, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, diagonal_letter, diagonal_positions, row, col + 1)

    if (row, col) in diagonal_positions:
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter, diagonal_positions, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diagonal_letter, diagonal_positions, row, col + 1):
                    return True
                grid[row][col] = ''

    return False

def solve_puzzle():
    grid = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]

    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

    for diagonal_letter in 'abcdefg':
        # Try each letter for the diagonal
        if solve(grid, diagonal_letter, diagonal_positions):
            result = '\n'.join(','.join(row) for row in grid)
            print(f"<<<\n{result}\n>>>")
            return

solve_puzzle()