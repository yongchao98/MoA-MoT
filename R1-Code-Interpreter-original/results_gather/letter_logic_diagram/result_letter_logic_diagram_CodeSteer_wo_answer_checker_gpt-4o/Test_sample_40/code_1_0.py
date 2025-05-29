def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter can be placed at grid[row][col]
    # Check row and column constraints
    if letter in grid[row]:
        return False
    if letter in (grid[i][col] for i in range(7)):
        return False
    # Check diagonal constraint
    if (row + col == 6) and (letter != diagonal_letter):
        return False
    return True

def solve(grid, row, col, diagonal_letter):
    if row == 7:
        # All rows are filled
        return True
    if col == 7:
        # Move to the next row
        return solve(grid, row + 1, 0, diagonal_letter)
    if grid[row][col] != ' ':
        # Skip pre-filled cells
        return solve(grid, row, col + 1, diagonal_letter)

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = ' '  # Backtrack

    return False

def solve_puzzle():
    grid = [
        [' ', 'e', ' ', 'g', ' ', 'c', 'b'],
        [' ', 'd', ' ', ' ', ' ', 'b', 'a'],
        [' ', ' ', ' ', 'c', ' ', ' ', ' '],
        [' ', ' ', 'c', 'b', ' ', 'e', 'd'],
        ['f', ' ', ' ', ' ', 'e', 'd', ' '],
        [' ', ' ', ' ', ' ', ' ', ' ', ' '],
        ['b', 'a', ' ', 'd', 'g', ' ', ' ']
    ]

    # Determine the letter for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in (' ', letter) for i in range(7)):
            diagonal_letter = letter
            break

    if solve(grid, 0, 0, diagonal_letter):
        result = '\n'.join(','.join(row) for row in grid)
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found.")

solve_puzzle()