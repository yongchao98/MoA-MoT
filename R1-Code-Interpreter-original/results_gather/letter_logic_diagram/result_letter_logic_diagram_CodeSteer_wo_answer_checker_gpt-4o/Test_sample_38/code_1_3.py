def solve_puzzle(grid):
    def is_valid(letter, row, col):
        # Check if the letter can be placed at grid[row][col]
        if letter in grid[row]:  # Check row
            return False
        if any(grid[r][col] == letter for r in range(7)):  # Check column
            return False
        return True

    def backtrack():
        # Find the next empty cell
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if is_valid(letter, row, col):
                            grid[row][col] = letter
                            print(f"Placing {letter} at ({row}, {col})")
                            if backtrack():
                                return True
                            print(f"Backtracking from ({row}, {col})")
                            grid[row][col] = ''  # Backtrack
                    return False
        return True

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(is_valid(letter, i, 6 - i) for i in range(7)):
            for i in range(7):
                grid[i][6 - i] = letter
            print(f"Trying letter '{letter}' for the minor diagonal")
            if backtrack():
                break
            for i in range(7):
                grid[i][6 - i] = ''

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid with given letters
initial_grid = [
    ['a', 'e', '', 'f', 'g', 'c', ''],
    ['', 'b', '', '', '', 'd', 'a'],
    ['b', 'f', 'g', 'c', 'd', 'a', 'e'],
    ['f', '', 'c', '', 'a', '', 'b'],
    ['g', 'c', 'd', 'a', 'e', 'b', 'f'],
    ['c', '', '', 'e', 'b', '', 'g'],
    ['d', '', '', '', 'f', 'g', 'c']
]

# Solve the puzzle
solution = solve_puzzle(initial_grid)
print(f"<<<\n{solution}\n>>>")