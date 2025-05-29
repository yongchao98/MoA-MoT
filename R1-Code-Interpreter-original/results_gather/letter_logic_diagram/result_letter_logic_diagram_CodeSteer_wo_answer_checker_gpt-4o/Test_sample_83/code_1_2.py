def solve_puzzle(grid):
    from collections import defaultdict

    # Helper function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    # Find the missing letters for each row and column
    rows_missing = [set('abcdefg') - set(row) for row in grid]
    cols_missing = [set('abcdefg') - set(grid[r][c] for r in range(7)) for c in range(7)]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(can_place(letter, r, 6-r) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
                rows_missing[r].discard(letter)
                cols_missing[6-r].discard(letter)

            # Forward checking and maintaining arc consistency
            def forward_check(row, col):
                if grid[row][col] != '':
                    return True
                for l in rows_missing[row] & cols_missing[col]:
                    grid[row][col] = l
                    rows_missing[row].remove(l)
                    cols_missing[col].remove(l)
                    if forward_check(row, col + 1 if col < 6 else 0):
                        return True
                    grid[row][col] = ''
                    rows_missing[row].add(l)
                    cols_missing[col].add(l)
                return False

            # Backtracking with constraint propagation
            def backtrack(row=0, col=0):
                if row == 7:
                    return True
                if col == 7:
                    return backtrack(row + 1, 0)
                if grid[row][col] != '':
                    return backtrack(row, col + 1)

                for l in rows_missing[row] & cols_missing[col]:
                    grid[row][col] = l
                    rows_missing[row].remove(l)
                    cols_missing[col].remove(l)
                    if forward_check(row, col + 1 if col < 6 else 0):
                        return True
                    grid[row][col] = ''
                    rows_missing[row].add(l)
                    cols_missing[col].add(l)

                return False

            if backtrack():
                return grid

    return None

# Initial grid setup
grid = [
    ['d', '', '', 'f', 'a', '', ''],
    ['b', '', 'f', '', '', 'e', 'd'],
    ['', 'f', 'a', 'c', 'e', 'd', 'b'],
    ['f', 'a', 'c', '', '', '', ''],
    ['a', '', 'e', '', '', 'g', 'f'],
    ['c', 'e', '', 'b', '', '', 'a'],
    ['e', 'd', '', 'g', 'f', 'a', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found.")