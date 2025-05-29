def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True  # All rows are filled

    if col == 7:
        return solve(grid, row + 1, 0)  # Move to the next row

    if grid[row][col] != '':
        return solve(grid, row, col + 1)  # Skip pre-filled cells

    # Try placing each letter in the current cell
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

def find_diagonal_letter(grid):
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.discard(grid[i][6-i])
    return possible_letters.pop()

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

    # Determine the letter for the minor diagonal
    diagonal_letter = find_diagonal_letter(grid)
    fill_diagonal(grid, diagonal_letter)

    # Solve the puzzle
    if solve(grid):
        result = '\n'.join([','.join(row) for row in grid])
        print(f"<<<\n{result}\n>>>")
    else:
        print("No solution found")

solve_puzzle()