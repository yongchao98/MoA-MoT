def solve_puzzle(grid):
    from collections import defaultdict

    def is_consistent(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    def select_unassigned_variable():
        # Use Minimum Remaining Values (MRV) heuristic
        min_options = float('inf')
        chosen_cell = None
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    options = len(rows_missing[r] & cols_missing[c])
                    if options < min_options:
                        min_options = options
                        chosen_cell = (r, c)
        return chosen_cell

    def order_domain_values(row, col):
        # Use Least Constraining Value (LCV) heuristic
        possible_letters = list(rows_missing[row] & cols_missing[col])
        constraints_count = defaultdict(int)
        for letter in possible_letters:
            for r in range(7):
                if letter in rows_missing[r]:
                    constraints_count[letter] += 1
            for c in range(7):
                if letter in cols_missing[c]:
                    constraints_count[letter] += 1
        return sorted(possible_letters, key=lambda l: constraints_count[l])

    def maintain_arc_consistency():
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    possible_letters = rows_missing[r] & cols_missing[c]
                    if not possible_letters:
                        return False
        return True

    def backtrack():
        if all(grid[r][c] != '' for r in range(7) for c in range(7)):
            return True

        row, col = select_unassigned_variable()
        for letter in order_domain_values(row, col):
            if is_consistent(letter, row, col):
                grid[row][col] = letter
                rows_missing[row].remove(letter)
                cols_missing[col].remove(letter)

                if maintain_arc_consistency() and backtrack():
                    return True

                # Backtrack
                grid[row][col] = ''
                rows_missing[row].add(letter)
                cols_missing[col].add(letter)

        return False

    # Find the missing letters for each row and column
    rows_missing = [set('abcdefg') - set(row) for row in grid]
    cols_missing = [set('abcdefg') - set(grid[r][c] for r in range(7)) for c in range(7)]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        if all(is_consistent(letter, r, 6-r) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
                rows_missing[r].discard(letter)
                cols_missing[6-r].discard(letter)

            # Try to fill the rest of the grid
            if backtrack():
                return grid

            # Reset the grid if unsuccessful
            for r in range(7):
                grid[r][6-r] = ''
                rows_missing[r].add(letter)
                cols_missing[6-r].add(letter)

    return None

# Initial grid
grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    print("<<<")
    for row in solution:
        print(','.join(row))
    print(">>>")
else:
    print("No solution found.")