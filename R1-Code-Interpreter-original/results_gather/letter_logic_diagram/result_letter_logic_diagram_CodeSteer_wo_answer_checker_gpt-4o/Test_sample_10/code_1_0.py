def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] is not None:
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = None  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['d', None, None, 'e', 'b', None, 'g'],
        [None, None, 'e', 'b', None, None, 'd'],
        [None, None, 'b', None, None, None, 'c'],
        [None, None, 'a', None, 'd', 'c', 'f'],
        [None, 'a', 'g', 'd', None, None, 'e'],
        [None, None, None, None, None, None, None],
        [None, None, None, 'f', None, None, 'a']
    ]

    # Fill the minor diagonal with a common letter
    for letter in 'abcdefg':
        if all(grid[i][6-i] is None or grid[i][6-i] == letter for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            break

    if solve(grid):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found")

solve_puzzle()