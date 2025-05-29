def solve_puzzle(grid):
    from itertools import permutations

    # Helper function to check if a letter can be placed at grid[row][col]
    def can_place(row, col, letter):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            minor_diagonal_letter = grid[i][6-i]
            break

    # If no pre-filled letter on the minor diagonal, choose one
    if not minor_diagonal_letter:
        for letter in 'abcdefg':
            if all(can_place(i, 6-i, letter) for i in range(7)):
                minor_diagonal_letter = letter
                break

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row):
        if row == 7:
            return True
        # Find the next empty cell in the current row
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if can_place(row, col, letter):
                        grid[row][col] = letter
                        if backtrack(row):
                            return True
                        grid[row][col] = ''  # Backtrack
                return False
        return backtrack(row + 1)

    # Start backtracking from the first row
    backtrack(0)

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid setup
grid = [
    ['', '', 'e', '', '', '', ''],
    ['g', 'e', 'f', '', '', '', 'c'],
    ['', '', '', 'a', '', '', ''],
    ['', '', '', '', '', '', 'e'],
    ['', 'a', 'd', 'c', '', '', ''],
    ['', 'd', 'c', 'g', '', '', 'b'],
    ['d', 'c', 'g', '', 'f', 'b', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")