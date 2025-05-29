def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve_puzzle(grid):
    def backtrack():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    for letter in 'abcdefg':
                        if is_valid(grid, i, j, letter):
                            grid[i][j] = letter
                            if backtrack():
                                return True
                            grid[i][j] = ''
                    return False
        return True

    # Fill the minor diagonal with the same letter
    diagonal_letter = 'a'
    for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        grid[i][j] = diagonal_letter

    if backtrack():
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "Error: No valid solution found"

# Initial grid
grid = [
    ['f', '', 'd', '', '', 'g', 'a'],
    ['', '', 'e', '', '', '', ''],
    ['', 'e', '', 'g', '', '', 'c'],
    ['', '', '', 'a', '', 'c', 'd'],
    ['', 'g', '', 'f', 'c', 'd', 'e'],
    ['', '', '', '', '', '', ''],
    ['', '', 'c', 'd', 'e', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)