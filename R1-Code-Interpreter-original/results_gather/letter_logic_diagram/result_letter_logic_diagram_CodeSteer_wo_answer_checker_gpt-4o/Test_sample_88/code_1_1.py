def solve_puzzle(grid):
    from collections import defaultdict

    # Initialize possible letters for each cell
    possible_letters = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Update possible letters based on initial grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                letter = grid[r][c]
                possible_letters[r][c] = {letter}
                # Remove this letter from possible letters in the same row and column
                for i in range(7):
                    if i != c:
                        possible_letters[r][i].discard(letter)
                    if i != r:
                        possible_letters[i][c].discard(letter)

    # Find the letter for the minor diagonal
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal with the chosen letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter
        possible_letters[i][6-i] = {diagonal_letter}

    # Function to check if the grid is complete
    def is_complete():
        return all(grid[r][c] != '' for r in range(7) for c in range(7))

    # Backtracking function with forward checking
    def backtrack():
        if is_complete():
            return True

        # Find the cell with the fewest possible letters
        min_options = 8
        cell = None
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '' and len(possible_letters[r][c]) < min_options:
                    min_options = len(possible_letters[r][c])
                    cell = (r, c)

        if cell is None:
            return False

        r, c = cell
        for letter in possible_letters[r][c]:
            # Try assigning the letter
            grid[r][c] = letter
            # Save current state of possible letters
            saved_state = [row[:] for row in possible_letters]
            # Update possible letters for the row and column
            for i in range(7):
                if i != c:
                    possible_letters[r][i].discard(letter)
                if i != r:
                    possible_letters[i][c].discard(letter)

            # Recurse
            if backtrack():
                return True

            # Restore state if failed
            grid[r][c] = ''
            possible_letters = saved_state

        return False

    # Start backtracking
    backtrack()

    # Format the output
    result = '\n'.join(','.join(row) for row in grid)
    return f"<<<\n{result}\n>>>"

# Initial grid setup
grid = [
    ['f', '', '', '', 'e', 'a', 'g'],
    ['c', 'd', 'b', '', 'a', 'g', 'f'],
    ['', 'b', '', 'a', 'g', 'f', 'c'],
    ['b', 'e', '', 'g', 'f', 'c', 'd'],
    ['', 'a', 'g', 'f', 'c', 'd', 'b'],
    ['', 'g', 'f', '', '', 'b', 'e'],
    ['g', '', '', '', '', 'e', 'a']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)