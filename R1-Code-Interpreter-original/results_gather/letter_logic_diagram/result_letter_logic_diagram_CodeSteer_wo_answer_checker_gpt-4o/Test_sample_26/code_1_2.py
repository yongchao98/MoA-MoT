def is_valid(grid, row, col, letter, diagonal_letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    # Check if the letter is valid for the diagonal
    if col == 6 - row and letter != diagonal_letter:
        return False
    return True

def solve(grid, row=0, col=0, diagonal_letter=''):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0, diagonal_letter)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1, diagonal_letter)  # Skip pre-filled cells

    # Try placing each letter in the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1, diagonal_letter):
                return True
            grid[row][col] = ''  # Backtrack

    return False

def fill_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def find_diagonal_letter(grid):
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.discard(grid[i][6-i])
    return possible_letters

def solve_puzzle():
    grid = [
        ['', 'g', 'a', '', 'd', '', ''],
        ['a', '', '', '', '', 'c', 'f'],
        ['a', '', 'd', 'b', 'c', '', 'g'],
        ['e', 'd', '', 'c', 'f', 'g', 'a'],
        ['', 'b', 'c', '', '', 'a', ''],
        ['', 'c', '', 'g', 'a', 'e', 'd'],
        ['', 'f', '', 'a', '', 'd', 'b']
    ]

    # Try each possible letter for the minor diagonal
    for diagonal_letter in find_diagonal_letter(grid):
        # Make a copy of the grid to try this diagonal letter
        test_grid = [row[:] for row in grid]
        fill_diagonal(test_grid, diagonal_letter)

        # Solve the puzzle
        if solve(test_grid, diagonal_letter=diagonal_letter):
            result = '\n'.join([','.join(row) for row in test_grid])
            print(f"<<<\n{result}\n>>>")
            return

    print("No solution found")

solve_puzzle()