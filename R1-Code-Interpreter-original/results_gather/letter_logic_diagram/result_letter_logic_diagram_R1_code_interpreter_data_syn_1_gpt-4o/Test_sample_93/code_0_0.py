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

    # Fill the minor diagonal with 'c'
    for i in range(7):
        grid[i][6-i] = 'c'

    # Try to fill the grid
    def fill_grid(row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return fill_grid(row + 1, 0)
        if grid[row][col] != '':
            return fill_grid(row, col + 1)

        for letter in 'abcdefg':
            if can_place(row, col, letter):
                grid[row][col] = letter
                if fill_grid(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    fill_grid()
    return grid

# Initial grid setup
initial_grid = [
    ['g', 'b', '', '', '', 'a', 'c'],
    ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
    ['', '', '', '', 'c', 'g', ''],
    ['f', '', '', '', '', '', 'e'],
    ['', 'a', 'c', '', 'b', 'e', 'f'],
    ['a', '', 'g', 'b', 'e', 'f', ''],
    ['c', 'g', 'b', '', '', 'd', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(initial_grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))