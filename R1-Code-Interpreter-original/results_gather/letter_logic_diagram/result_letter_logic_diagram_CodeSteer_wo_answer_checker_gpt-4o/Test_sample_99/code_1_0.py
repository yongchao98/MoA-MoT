def solve_puzzle(grid):
    from collections import deque
    from copy import deepcopy

    # Determine possible letters for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_diagonal_letters = set('abcdefg')
    for i, j in diagonal_indices:
        if grid[i][j] != '':
            possible_diagonal_letters.intersection_update(grid[i][j])

    # Function to initialize possible values for each cell
    def initialize_possibilities():
        possibilities = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
        for r in range(7):
            for c in range(7):
                if grid[r][c] != '':
                    possibilities[r][c] = set(grid[r][c])
        return possibilities

    # Function to update possibilities after placing a letter
    def update_possibilities(possibilities, row, col, letter):
        for c in range(7):
            possibilities[row][c].discard(letter)
        for r in range(7):
            possibilities[r][col].discard(letter)
        for i, j in diagonal_indices:
            possibilities[i][j].discard(letter)

    # AC-3 algorithm to maintain arc consistency
    def ac3(possibilities):
        queue = deque([(r, c) for r in range(7) for c in range(7) if len(possibilities[r][c]) == 1])
        while queue:
            row, col = queue.popleft()
            if len(possibilities[row][col]) == 1:
                letter = next(iter(possibilities[row][col]))
                for c in range(7):
                    if c != col and letter in possibilities[row][c]:
                        possibilities[row][c].discard(letter)
                        if len(possibilities[row][c]) == 1:
                            queue.append((row, c))
                for r in range(7):
                    if r != row and letter in possibilities[r][col]:
                        possibilities[r][col].discard(letter)
                        if len(possibilities[r][col]) == 1:
                            queue.append((r, col))
                for i, j in diagonal_indices:
                    if (i, j) != (row, col) and letter in possibilities[i][j]:
                        possibilities[i][j].discard(letter)
                        if len(possibilities[i][j]) == 1:
                            queue.append((i, j))

    # Function to find the cell with the fewest possibilities
    def select_unassigned_variable(possibilities):
        min_possibilities = float('inf')
        best_cell = None
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '' and 0 < len(possibilities[r][c]) < min_possibilities:
                    min_possibilities = len(possibilities[r][c])
                    best_cell = (r, c)
        return best_cell

    # Backtracking function with constraint propagation
    def backtrack(grid, possibilities):
        ac3(possibilities)
        cell = select_unassigned_variable(possibilities)
        if not cell:
            return True
        row, col = cell

        for letter in possibilities[row][col]:
            grid[row][col] = letter
            new_possibilities = deepcopy(possibilities)
            update_possibilities(new_possibilities, row, col, letter)
            if backtrack(grid, new_possibilities):
                return True
            grid[row][col] = ''

        return False

    # Try each possible letter for the diagonal
    for diagonal_letter in possible_diagonal_letters:
        # Fill the diagonal with the current letter
        for i, j in diagonal_indices:
            grid[i][j] = diagonal_letter

        # Initialize possibilities
        possibilities = initialize_possibilities()

        # Start backtracking from the first cell
        if backtrack(grid, possibilities):
            # Format the output
            result = '\n'.join(','.join(row) for row in grid)
            return f"<<<\n{result}\n>>>"

    return "No solution found"

# Initial grid with empty cells as ''
grid = [
    ['', '', 'b', '', '', 'd', 'a'],
    ['', '', 'e', '', 'd', '', 'g'],
    ['', 'e', 'c', 'd', '', '', ''],
    ['', 'c', '', '', '', '', ''],
    ['', 'd', 'a', '', 'f', '', 'e'],
    ['', '', 'g', '', '', '', ''],
    ['', '', 'f', '', 'e', 'c', 'd']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)