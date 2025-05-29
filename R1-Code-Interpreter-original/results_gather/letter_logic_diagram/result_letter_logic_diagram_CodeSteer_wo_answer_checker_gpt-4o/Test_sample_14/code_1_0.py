def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter is already in the current row or column
        if letter in rows[row] or letter in cols[col]:
            return False
        # Check the minor diagonal constraint
        if row + col == 6 and letter != minor_diagonal_letter:
            return False
        return True

    def solve(grid, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return solve(grid, row + 1, 0)
        if grid[row][col] != '':
            return solve(grid, row, col + 1)

        for letter in possible_values[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                rows[row].add(letter)
                cols[col].add(letter)
                if row + col == 6:
                    diagonal.add(letter)

                if solve(grid, row, col + 1):
                    return True

                grid[row][col] = ''
                rows[row].remove(letter)
                cols[col].remove(letter)
                if row + col == 6:
                    diagonal.remove(letter)
        return False

    # Convert input grid to a list of lists
    grid = [list(row.split(',')) for row in grid.strip().split('\n')]

    # Initialize sets for rows, columns, and diagonal
    rows = [set() for _ in range(7)]
    cols = [set() for _ in range(7)]
    diagonal = set()

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6 - i] != '':
            minor_diagonal_letter = grid[i][6 - i]
            diagonal.add(minor_diagonal_letter)
            break
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'  # Default to 'a' if no pre-filled letter

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        if grid[i][6 - i] == '':
            grid[i][6 - i] = minor_diagonal_letter
        rows[i].add(grid[i][6 - i])
        cols[6 - i].add(grid[i][6 - i])

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                letter = grid[r][c]
                rows[r].add(letter)
                cols[c].add(letter)
                possible_values[r][c] = {letter}

    # Solve the puzzle
    if solve(grid):
        # Print the solution
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Input grid
puzzle = """
,f,,c,,,
,d,,e,,,
,,,,,b,
c,e,g,,b,f,
,,,,f,d,
,a,b,,,c,
,b,f,d,,,g
"""

solve_puzzle(puzzle)