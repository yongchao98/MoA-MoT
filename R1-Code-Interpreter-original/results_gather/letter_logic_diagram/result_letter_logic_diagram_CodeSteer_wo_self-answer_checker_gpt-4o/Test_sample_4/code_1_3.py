def is_valid(grid, row, col, letter, row_sets, col_sets):
    # Check if the letter is already in the row or column
    if letter in row_sets[row] or letter in col_sets[col]:
        return False
    return True

def solve(grid, diagonal_letter, row_sets, col_sets, row=0, col=0):
    if row == 7:
        return True  # Puzzle solved

    if col == 7:
        return solve(grid, diagonal_letter, row_sets, col_sets, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, diagonal_letter, row_sets, col_sets, row, col + 1)  # Skip pre-filled cells

    # Determine if this cell is on the minor diagonal
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter, row_sets, col_sets):
            grid[row][col] = diagonal_letter
            row_sets[row].add(diagonal_letter)
            col_sets[col].add(diagonal_letter)
            if solve(grid, diagonal_letter, row_sets, col_sets, row, col + 1):
                return True
            grid[row][col] = ''  # Backtrack
            row_sets[row].remove(diagonal_letter)
            col_sets[col].remove(diagonal_letter)
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, row_sets, col_sets):
                grid[row][col] = letter
                row_sets[row].add(letter)
                col_sets[col].add(letter)
                if solve(grid, diagonal_letter, row_sets, col_sets, row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack
                row_sets[row].remove(letter)
                col_sets[col].remove(letter)

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

    # Initialize sets to track used letters in each row and column
    row_sets = [set(row) for row in grid]
    col_sets = [set(grid[i][j] for i in range(7) if grid[i][j]) for j in range(7)]

    # Solve the puzzle
    if solve(grid, diagonal_letter, row_sets, col_sets):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()