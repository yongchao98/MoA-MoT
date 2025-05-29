def solve_puzzle(grid):
    # Initialize sets to track available letters for each row and column
    rows = [set('abcdefg') for _ in range(7)]
    cols = [set('abcdefg') for _ in range(7)]

    # Fill the sets with the pre-filled letters
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                rows[r].discard(grid[r][c])
                cols[c].discard(grid[r][c])

    def is_valid(letter, row, col):
        return letter in rows[row] and letter in cols[col]

    def backtrack():
        # Find the next empty cell with the fewest options
        min_options = 8
        next_cell = None
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    options = rows[r] & cols[c]
                    if len(options) < min_options:
                        min_options = len(options)
                        next_cell = (r, c)
        
        if not next_cell:
            return True  # All cells are filled

        row, col = next_cell
        for letter in rows[row] & cols[col]:
            grid[row][col] = letter
            rows[row].remove(letter)
            cols[col].remove(letter)

            if backtrack():
                return True

            # Backtrack
            grid[row][col] = ''
            rows[row].add(letter)
            cols[col].add(letter)

        return False

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(is_valid(letter, i, 6 - i) for i in range(7)):
            for i in range(7):
                grid[i][6 - i] = letter
                rows[i].discard(letter)
                cols[6 - i].discard(letter)

            if backtrack():
                break

            # Backtrack the minor diagonal
            for i in range(7):
                grid[i][6 - i] = ''
                rows[i].add(letter)
                cols[6 - i].add(letter)

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid setup
grid = [
    ['c', '', '', '', 'b', '', ''],
    ['', 'f', 'e', 'b', '', 'g', ''],
    ['', '', 'b', '', '', '', ''],
    ['', 'b', 'a', '', 'c', 'd', 'f'],
    ['b', 'a', '', '', '', '', ''],
    ['', 'g', '', '', 'f', 'e', ''],
    ['', 'c', '', 'f', '', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")