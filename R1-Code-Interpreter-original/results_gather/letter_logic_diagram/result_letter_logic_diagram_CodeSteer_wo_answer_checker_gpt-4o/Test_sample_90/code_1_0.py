def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['c', 'b', '', 'a', 'f', 'g', ''],
        ['b', 'e', 'a', 'f', '', '', ''],
        ['', '', '', 'g', '', '', ''],
        ['', '', '', '', '', 'b', ''],
        ['', 'g', '', 'c', 'b', '', ''],
        ['g', 'd', '', '', '', '', 'f'],
        ['', '', 'b', '', 'a', '', '']
    ]

    # Determine the minor diagonal letter
    minor_diagonal_letter = 'd'  # Based on analysis of pre-filled diagonal positions
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter

    if solve(grid):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()