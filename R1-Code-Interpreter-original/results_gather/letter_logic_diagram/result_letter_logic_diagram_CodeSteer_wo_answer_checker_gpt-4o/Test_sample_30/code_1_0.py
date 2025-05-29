def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Solved

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    # Try placing each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def fill_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def solve_puzzle():
    grid = [
        ['', '', 'a', '', 'c', 'f', 'b'],
        ['g', 'a', '', '', '', 'b', 'e'],
        ['', '', 'c', 'f', 'b', 'e', 'g'],
        ['d', 'c', 'f', '', 'e', 'g', 'a'],
        ['', 'f', 'b', 'e', 'g', '', ''],
        ['f', 'b', 'e', 'g', 'a', 'd', ''],
        ['b', 'e', '', 'a', 'd', '', 'f']
    ]

    # Try each letter for the diagonal
    for letter in 'abcdefg':
        fill_diagonal(grid, letter)
        if solve(grid):
            break

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    print(f"<<<\n{result}\n>>>")

solve_puzzle()