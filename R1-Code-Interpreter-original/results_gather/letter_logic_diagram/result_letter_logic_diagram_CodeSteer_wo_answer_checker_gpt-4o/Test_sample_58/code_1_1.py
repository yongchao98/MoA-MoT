def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def forward_check(grid, possibilities, row, col, letter):
        # Update possibilities after placing a letter
        for i in range(7):
            if letter in possibilities[row][i]:
                possibilities[row][i].remove(letter)
            if letter in possibilities[i][col]:
                possibilities[i][col].remove(letter)

    def fill_grid(grid, possibilities, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return fill_grid(grid, possibilities, row + 1, 0)
        if grid[row][col] != '':
            return fill_grid(grid, possibilities, row, col + 1)

        for letter in possibilities[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                new_possibilities = [row[:] for row in possibilities]
                forward_check(grid, new_possibilities, row, col, letter)
                if fill_grid(grid, new_possibilities, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Initialize possibilities
    possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                possibilities[r][c] = set(grid[r][c])
                forward_check(grid, possibilities, r, c, grid[r][c])

    # Fill the minor diagonal with the same letter
    for letter in 'abcdefg':
        temp_grid = [row[:] for row in grid]
        temp_possibilities = [row[:] for row in possibilities]
        for i in range(7):
            temp_grid[i][6-i] = letter
            temp_possibilities[i][6-i] = set(letter)
        if fill_grid(temp_grid, temp_possibilities):
            return temp_grid

    return None

# Initial grid setup
grid = [
    ['d', '', '', '', '', 'f', ''],
    ['', 'b', '', '', 'f', 'a', 'd'],
    ['b', '', '', 'f', '', '', 'd'],
    ['', '', 'f', '', '', '', ''],
    ['c', '', '', '', 'e', '', 'g'],
    ['', 'a', 'd', '', '', '', 'c'],
    ['', '', 'e', '', '', 'c', 'f']
]

solution = solve_puzzle(grid)
if solution:
    print("<<<")
    for row in solution:
        print(','.join(row))
    print(">>>")
else:
    print("No solution found.")