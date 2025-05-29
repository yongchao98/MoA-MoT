def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Solved

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
        grid[row][col] = ''  # Backtrack

    return False

# Initial grid setup
grid = [
    ['a', 'e', 'f', 'g', '', 'c', 'd'],
    ['e', 'f', 'g', 'b', '', 'd', 'a'],
    ['', 'g', 'b', '', '', '', 'e'],
    ['', 'b', '', '', 'a', '', 'f'],
    ['', 'c', 'd', '', 'e', '', 'g'],
    ['c', 'd', 'a', 'e', 'f', 'g', 'b'],
    ['d', 'a', 'e', 'f', '', 'b', '']
]

# Try each letter for the minor diagonal
for minor_diagonal_letter in 'abcdefg':
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
        break
else:
    print("No solution found")