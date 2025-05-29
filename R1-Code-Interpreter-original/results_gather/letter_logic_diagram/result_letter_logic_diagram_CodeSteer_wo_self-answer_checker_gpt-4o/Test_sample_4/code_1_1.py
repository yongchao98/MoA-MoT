def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Puzzle solved

    if col == 7:
        return solve(grid, diagonal_letter, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, diagonal_letter, row, col + 1)  # Skip pre-filled cells

    # Determine if this cell is on the minor diagonal
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['b', '', '', '', 'g', '', 'e'],
        ['', '', 'a', 'g', '', '', ''],
        ['f', 'a', '', '', '', 'b', 'd'],
        ['', 'g', '', '', '', '', ''],
        ['g', '', '', 'b', '', '', ''],
        ['', 'e', 'b', 'd', 'f', 'a', ''],
        ['', 'b', '', '', 'a', '', '']
    ]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_positions if grid[i][j])
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Solve the puzzle
    if solve(grid, diagonal_letter):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()