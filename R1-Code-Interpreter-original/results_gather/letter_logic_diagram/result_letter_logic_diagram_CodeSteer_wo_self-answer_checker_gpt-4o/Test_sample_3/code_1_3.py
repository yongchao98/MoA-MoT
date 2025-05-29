def solve_puzzle(grid):
    def is_valid(row, col, char):
        return char in rows[row] and char in cols[col]

    def place_letter(row, col, char):
        grid[row][col] = char
        rows[row].remove(char)
        cols[col].remove(char)

    def remove_letter(row, col, char):
        grid[row][col] = ''
        rows[row].add(char)
        cols[col].add(char)

    def solve(row, col):
        if row == 7:
            return True
        if col == 7:
            return solve(row + 1, 0)
        if grid[row][col] != '':
            return solve(row, col + 1)

        for char in 'abcdefg':
            if is_valid(row, col, char):
                place_letter(row, col, char)
                if solve(row, col + 1):
                    return True
                remove_letter(row, col, char)
        return False

    # Initialize sets for available letters in each row and column
    rows = [set('abcdefg') for _ in range(7)]
    cols = [set('abcdefg') for _ in range(7)]

    # Fill the minor diagonal with the same letter
    minor_diagonal_letter = 'g'  # From the given grid, 'g' is on the diagonal
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    # Remove pre-filled letters from available sets
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                if grid[r][c] in rows[r]:
                    rows[r].remove(grid[r][c])
                if grid[r][c] in cols[c]:
                    cols[c].remove(grid[r][c])

    # Solve the puzzle
    solve(0, 0)

    # Print the solution
    for row in grid:
        print(','.join(row))

# Initial grid with pre-filled values
grid = [
    ['g', '', 'c', 'a', 'd', '', 'e'],
    ['', 'c', '', '', '', '', ''],
    ['', '', '', '', '', 'g', ''],
    ['', '', 'b', '', 'g', 'f', ''],
    ['d', '', '', '', 'f', '', 'a'],
    ['b', '', 'g', 'f', '', 'a', ''],
    ['', '', '', '', 'a', 'd', '']
]

solve_puzzle(grid)