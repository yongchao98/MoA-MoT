def solve_puzzle(grid):
    from collections import defaultdict

    # Initialize possible candidates for each cell
    candidates = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Update candidates based on the initial grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                candidates[r][c] = {grid[r][c]}
                # Remove this letter from the row and column candidates
                for i in range(7):
                    if i != c:
                        candidates[r][i].discard(grid[r][c])
                    if i != r:
                        candidates[i][c].discard(grid[r][c])

    # Function to update candidates after placing a letter
    def update_candidates(r, c, letter):
        for i in range(7):
            if i != c:
                candidates[r][i].discard(letter)
            if i != r:
                candidates[i][c].discard(letter)

    # Function to select the next cell to fill using MRV heuristic
    def select_next_cell():
        min_candidates = 8  # More than the max possible (7)
        next_cell = None
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '' and 0 < len(candidates[r][c]) < min_candidates:
                    min_candidates = len(candidates[r][c])
                    next_cell = (r, c)
        return next_cell

    # Try to fill the minor diagonal with a single letter
    for letter in 'abcdefg':
        if all(letter in candidates[i][6-i] for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
                update_candidates(i, 6-i, letter)
            break

    # Backtracking function with constraint propagation
    def backtrack():
        cell = select_next_cell()
        if not cell:
            return True
        r, c = cell
        for letter in list(candidates[r][c]):
            grid[r][c] = letter
            original_candidates = [candidates[r][i].copy() for i in range(7)] + [candidates[i][c].copy() for i in range(7)]
            update_candidates(r, c, letter)
            if backtrack():
                return True
            # Restore candidates if backtracking
            for i in range(7):
                if i != c:
                    candidates[r][i] = original_candidates[i]
                if i != r:
                    candidates[i][c] = original_candidates[7 + i]
            grid[r][c] = ''
        return False

    backtrack()

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid setup
grid = [
    ['b', '', '', 'g', '', 'c', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', '', 'a', 'c', '', 'b', 'f'],
    ['g', '', '', 'e', '', '', 'd'],
    ['', '', 'e', 'b', '', '', ''],
    ['', 'e', '', 'f', 'd', 'g', 'a'],
    ['e', 'b', 'f', 'd', 'g', '', 'c']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")