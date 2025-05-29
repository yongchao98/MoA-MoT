def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter can be placed at grid[row][col]
    # Check row and column constraints
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check diagonal constraint
    if col == 6 - row and letter != diagonal_letter:
        return False
    return True

def solve(grid, row, col, diagonal_letter):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1, diagonal_letter)  # Skip pre-filled cells

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['','e','','','','g',''],
        ['e','c','a','b','','d','f'],
        ['','b','g','','','f','e'],
        ['','b','g','','','e','c'],
        ['b','','d','f','','c',''],
        ['g','d','','','c','a','b'],
        ['d','','e','','a','b','']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letter = grid[i][6-i]
            break

    if not diagonal_letter:
        # Choose a letter that can fit without violating constraints
        for letter in 'abcdefg':
            if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
                diagonal_letter = letter
                break

    # Pre-fill the diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    if solve(grid, 0, 0, diagonal_letter):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found")

solve_puzzle()