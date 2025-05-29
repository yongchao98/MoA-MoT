def solve_puzzle(grid):
    import copy

    # Helper function to find missing letters in a list
    def find_missing_letters(lst):
        return [chr(i) for i in range(ord('a'), ord('g') + 1) if chr(i) not in lst]

    # Convert grid to a list of lists
    grid = [row.split(',') for row in grid.strip().split('\n')]

    # Determine the letter for the minor diagonal
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for r, c in diagonal_positions:
        if grid[r][c] != '':
            diagonal_letter = grid[r][c]
            break

    if not diagonal_letter:
        for letter in 'abcdefg':
            if all(grid[r][c] == '' or grid[r][c] == letter for r, c in diagonal_positions):
                diagonal_letter = letter
                break

    # Fill the diagonal with the chosen letter
    for r, c in diagonal_positions:
        grid[r][c] = diagonal_letter

    # Helper function to check if placing a letter is valid
    def is_valid(r, c, letter):
        # Check row
        if letter in grid[r]:
            return False
        # Check column
        if letter in [grid[i][c] for i in range(7)]:
            return False
        return True

    # Backtracking function
    def backtrack(r, c):
        if r == 7:
            return True
        if c == 7:
            return backtrack(r + 1, 0)
        if grid[r][c] != '':
            return backtrack(r, c + 1)

        for letter in 'abcdefg':
            if is_valid(r, c, letter):
                grid[r][c] = letter
                if backtrack(r, c + 1):
                    return True
                grid[r][c] = ''

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return f"<<<\n{result}\n>>>"

# Given puzzle
puzzle = """
f,d,,,,c,g
d,e,,a,c,,f
e,b,,,g,f,
,,c,g,f,,e
a,c,g,f,d,e,b
,g,,d,e,,a
g,f,,e,b,a,c
"""

# Solve the puzzle
print(solve_puzzle(puzzle))