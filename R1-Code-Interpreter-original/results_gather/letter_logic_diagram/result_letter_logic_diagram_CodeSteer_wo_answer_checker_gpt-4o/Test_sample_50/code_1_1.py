def solve_puzzle(grid):
    def is_valid(letter, row, col):
        # Check if the letter can be placed in the current row and column
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if is_valid(letter, row, col):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(is_valid(letter, i, 6 - i) for i in range(7)):
            for i in range(7):
                grid[i][6 - i] = letter
            if backtrack(0, 0):
                break
            for i in range(7):
                grid[i][6 - i] = ''

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid setup
grid = [
    ['c', '', '', '', 'b', '', ''],
    ['', 'f', 'e', 'b', '', 'g', ''],
    ['', '', 'b', '', '', '', ''],
    ['', 'b', 'a', '', 'c', 'd', 'f'],
    ['b', 'a', '', '', '', '', ''],
    ['', 'g', '', '', 'f', 'e', ''],
    ['', 'c', '', 'f', '', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")