def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def forward_checking(grid, row, col, letter):
        # Remove the letter from the domains of unassigned cells in the same row and column
        for i in range(7):
            if grid[row][i] == '':
                domains[row][i].discard(letter)
            if grid[i][col] == '':
                domains[i][col].discard(letter)

    def restore_domains(grid, row, col, letter):
        # Restore the letter to the domains of unassigned cells in the same row and column
        for i in range(7):
            if grid[row][i] == '':
                domains[row][i].add(letter)
            if grid[i][col] == '':
                domains[i][col].add(letter)

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in list(domains[row][col]):
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                forward_checking(grid, row, col, letter)
                if solve(grid, row, col + 1):
                    return True
                grid[row][col] = ''
                restore_domains(grid, row, col, letter)
        return False

    # Initialize the grid and domains
    grid = [list(row) for row in grid]
    domains = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Update domains based on the initial grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] != '':
                letter = grid[row][col]
                domains[row][col] = {letter}
                forward_checking(grid, row, col, letter)

    # Solve the puzzle
    if solve(grid):
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Given puzzle
puzzle = [
    ['', 'b', '', '', '', '', ''],
    ['b', 'a', 'g', 'd', '', '', 'f'],
    ['', '', 'd', '', '', '', 'b'],
    ['g', '', '', 'e', '', '', ''],
    ['', '', '', '', 'b', '', ''],
    ['c', '', 'f', '', 'a', '', 'd'],
    ['e', 'f', '', 'a', 'g', '', '']
]

# Solve and print the solution
solution = solve_puzzle(puzzle)
print("<<<")
print(solution)
print(">>>")