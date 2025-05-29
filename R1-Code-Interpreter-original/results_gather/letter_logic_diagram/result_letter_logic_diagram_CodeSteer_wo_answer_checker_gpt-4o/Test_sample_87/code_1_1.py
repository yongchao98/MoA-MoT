def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row, col, minor_diagonal_letter):
    if row == 7:
        return True  # Successfully filled the grid

    if col == 7:
        return solve(grid, row + 1, 0, minor_diagonal_letter)

    if grid[row][col] != '':
        return solve(grid, row, col + 1, minor_diagonal_letter)

    for letter in 'abcdefg':
        if (col == 6 - row and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue

        grid[row][col] = letter
        if solve(grid, row, col + 1, minor_diagonal_letter):
            return True
        grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle(puzzle):
    grid = [row.split(',') for row in puzzle.strip().split('\n')]
    letters = set('abcdefg')

    for minor_diagonal_letter in letters:
        # Try each letter as the minor diagonal letter
        temp_grid = [row[:] for row in grid]  # Copy the grid
        for i in range(7):
            if temp_grid[i][6-i] == '' or temp_grid[i][6-i] == minor_diagonal_letter:
                temp_grid[i][6-i] = minor_diagonal_letter
            else:
                break
        else:
            if solve(temp_grid, 0, 0, minor_diagonal_letter):
                result = '\n'.join(','.join(row) for row in temp_grid)
                return f"<<<\n{result}\n>>>"

    return "No valid configuration found."

# The given puzzle
puzzle = """
e,a,c,,,,
,,g,f,,d,e
c,g,f,b,d,,
,,b,d,e,,c
,,,,,,
,d,e,,c,g,f
d,,,,g,f,
"""

# Solve the puzzle
print(solve_puzzle(puzzle))