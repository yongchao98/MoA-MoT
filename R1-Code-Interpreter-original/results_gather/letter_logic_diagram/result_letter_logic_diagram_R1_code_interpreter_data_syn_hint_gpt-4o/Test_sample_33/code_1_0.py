def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter):
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        for letter in 'abcdefg':
            if (col == 6 - row and letter != diagonal_letter) or not is_valid(grid, row, col, letter):
                continue
            grid[row][col] = letter
            if backtrack(row, col + 1):
                return True
            grid[row][col] = ''
        
        return False

    return backtrack(0, 0)

# Initial grid
grid = [
    ['a', 'b', '', 'g', 'e', 'c', 'd'],
    ['', 'f', 'g', 'e', '', 'd', ''],
    ['f', 'g', 'e', '', '', 'a', ''],
    ['', 'e', 'c', 'd', '', '', 'f'],
    ['e', '', 'd', '', '', 'f', ''],
    ['c', '', 'a', 'b', '', 'g', ''],
    ['', '', 'b', '', 'g', 'e', 'c']
]

# Determine the letter for the minor diagonal
diagonal_letter = None
for letter in 'abcdefg':
    if all(grid[i][6-i] in ('', letter) for i in range(7)):
        diagonal_letter = letter
        break

# Solve the grid
solve(grid, diagonal_letter)

# Print the filled grid
for row in grid:
    print(','.join(row))