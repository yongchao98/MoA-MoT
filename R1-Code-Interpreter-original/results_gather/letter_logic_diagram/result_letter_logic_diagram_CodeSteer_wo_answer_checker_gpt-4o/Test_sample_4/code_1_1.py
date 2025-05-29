def solve_puzzle(grid):
    from collections import defaultdict

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Fill the minor diagonal with the same letter
    minor_diagonal = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letter = None
    for i, j in minor_diagonal:
        if grid[i][j] != '':
            diagonal_letter = grid[i][j]
            break
    if not diagonal_letter:
        # If no pre-filled letter, choose one that is not over-represented
        letter_counts = defaultdict(int)
        for row in grid:
            for letter in row:
                if letter:
                    letter_counts[letter] += 1
        diagonal_letter = min(letter_counts, key=letter_counts.get)

    for i, j in minor_diagonal:
        grid[i][j] = diagonal_letter
        possible_values[i][j] = {diagonal_letter}

    # Update possible values based on the current grid
    def update_possible_values():
        for i in range(7):
            for j in range(7):
                if grid[i][j] != '':
                    possible_values[i][j] = {grid[i][j]}
                else:
                    row_values = {grid[i][k] for k in range(7) if grid[i][k] != ''}
                    col_values = {grid[k][j] for k in range(7) if grid[k][j] != ''}
                    possible_values[i][j] -= row_values | col_values

    update_possible_values()

    # Backtracking function
    def backtrack():
        # Find the cell with the fewest possibilities
        min_possibilities = 8
        cell = None
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and len(possible_values[i][j]) < min_possibilities:
                    min_possibilities = len(possible_values[i][j])
                    cell = (i, j)
        
        if not cell:
            return True  # All cells are filled

        i, j = cell
        for value in possible_values[i][j]:
            grid[i][j] = value
            original_possibilities = [row[:] for row in possible_values]
            update_possible_values()
            if backtrack():
                return True
            grid[i][j] = ''
            possible_values[:] = original_possibilities  # Restore possibilities

        return False

    if backtrack():
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Initial grid setup
grid = [
    ['b', '', '', '', 'g', '', 'e'],
    ['', '', 'a', 'g', '', '', ''],
    ['f', 'a', '', '', '', 'b', 'd'],
    ['', 'g', '', '', '', '', ''],
    ['g', '', '', 'b', '', '', ''],
    ['', 'e', 'b', 'd', 'f', 'a', ''],
    ['', 'b', '', '', 'a', '', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")