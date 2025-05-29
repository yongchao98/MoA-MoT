def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def update_possibilities(possibilities, row, col, letter):
        # Remove the letter from the possibilities of the same row and column
        for i in range(7):
            possibilities[row][i].discard(letter)
            possibilities[i][col].discard(letter)

    def solve(grid, possibilities, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, possibilities, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, possibilities, row, col + 1)

        for letter in possibilities[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                new_possibilities = [row.copy() for row in possibilities]
                update_possibilities(new_possibilities, row, col, letter)
                if solve(grid, new_possibilities, row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    # Initialize possibilities
    possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                possibilities[r][c] = {grid[r][c]}
                update_possibilities(possibilities, r, c, grid[r][c])

    # Fill the minor diagonal with the same letter
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            if solve(grid, possibilities):
                return grid
            for i in range(7):
                grid[i][6-i] = ''

    return None

# Initial grid
grid = [
    ['', 'c', 'g', 'b', 'a', 'e', 'f'],
    ['', '', '', 'a', 'e', 'f', 'd'],
    ['g', 'b', '', 'e', '', '', 'c'],
    ['b', 'a', 'e', '', 'd', 'c', ''],
    ['a', 'e', 'f', 'd', 'c', '', 'b'],
    ['e', 'f', 'd', '', 'g', 'b', ''],
    ['f', '', '', 'g', '', 'a', 'e']
]

solution = solve_puzzle(grid)
if solution:
    print("<<<")
    for row in solution:
        print(','.join(row))
    print(">>>")
else:
    print("No solution found.")