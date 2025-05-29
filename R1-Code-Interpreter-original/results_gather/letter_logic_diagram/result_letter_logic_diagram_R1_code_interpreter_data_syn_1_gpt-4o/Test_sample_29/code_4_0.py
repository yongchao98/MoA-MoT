def solve_puzzle(grid):
    from collections import Counter

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the missing letters for each row and column
    all_letters = set('abcdefg')
    missing_in_rows = [all_letters - set(row) for row in grid]
    missing_in_cols = [all_letters - set(col) for col in zip(*grid)]

    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(i, 6-i) for i in range(7)]
    minor_diagonal_letters = [grid[i][6-i] for i in range(7) if grid[i][6-i]]
    if minor_diagonal_letters:
        minor_diagonal_letter = Counter(minor_diagonal_letters).most_common(1)[0][0]
    else:
        # If no letters are pre-filled on the diagonal, choose one that can fit
        possible_letters = all_letters - set(minor_diagonal_letters)
        minor_diagonal_letter = possible_letters.pop()

    # Function to check if placing a letter is valid
    def is_valid(letter, row, col):
        return letter in missing_in_rows[row] and letter in missing_in_cols[col]

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col]:
            return backtrack(row, col + 1)

        if (row, col) in minor_diagonal_indices:
            if is_valid(minor_diagonal_letter, row, col):
                grid[row][col] = minor_diagonal_letter
                missing_in_rows[row].remove(minor_diagonal_letter)
                missing_in_cols[col].remove(minor_diagonal_letter)
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
                missing_in_rows[row].add(minor_diagonal_letter)
                missing_in_cols[col].add(minor_diagonal_letter)
        else:
            for letter in all_letters:
                if is_valid(letter, row, col):
                    grid[row][col] = letter
                    missing_in_rows[row].remove(letter)
                    missing_in_cols[col].remove(letter)
                    if backtrack(row, col + 1):
                        return True
                    grid[row][col] = ''
                    missing_in_rows[row].add(letter)
                    missing_in_cols[col].add(letter)

        return False

    # Start backtracking from the first cell
    if backtrack(0, 0):
        return '\n'.join([','.join(row) for row in grid])
    else:
        return "No solution found."

# Given puzzle
puzzle = """
a,e,f,g,,c,d
e,f,g,b,,d,a
,g,b,,,,e
,b,,,a,,f
,c,d,,e,,g
c,d,a,e,f,g,b
d,a,e,f,,b,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(solution)