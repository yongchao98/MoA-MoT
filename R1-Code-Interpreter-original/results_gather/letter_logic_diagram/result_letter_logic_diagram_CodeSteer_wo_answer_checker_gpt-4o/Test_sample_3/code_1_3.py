def solve_puzzle(grid):
    from collections import deque

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]

    # Update possible values based on the initial grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] != '':
                letter = grid[r][c]
                possible_values[r][c] = {letter}
                for i in range(7):
                    possible_values[r][i].discard(letter)
                    possible_values[i][c].discard(letter)

    # Function to enforce arc consistency (AC-3)
    def ac3():
        queue = deque((r, c) for r in range(7) for c in range(7) if grid[r][c] == '')
        while queue:
            r, c = queue.popleft()
            if len(possible_values[r][c]) == 1:
                letter = next(iter(possible_values[r][c]))
                for i in range(7):
                    if i != c and letter in possible_values[r][i]:
                        possible_values[r][i].discard(letter)
                        if len(possible_values[r][i]) == 1:
                            queue.append((r, i))
                    if i != r and letter in possible_values[i][c]:
                        possible_values[i][c].discard(letter)
                        if len(possible_values[i][c]) == 1:
                            queue.append((i, c))

    # Function to select the next cell to fill using MRV and Degree Heuristic
    def select_unassigned_cell():
        min_possibilities = 8
        best_cell = None
        max_degree = -1
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    num_possibilities = len(possible_values[r][c])
                    if num_possibilities < min_possibilities:
                        min_possibilities = num_possibilities
                        best_cell = (r, c)
                        max_degree = sum(1 for i in range(7) if grid[r][i] == '' or grid[i][c] == '')
                    elif num_possibilities == min_possibilities:
                        degree = sum(1 for i in range(7) if grid[r][i] == '' or grid[i][c] == '')
                        if degree > max_degree:
                            best_cell = (r, c)
                            max_degree = degree
        return best_cell

    # Backtracking function with forward checking
    def backtrack():
        ac3()  # Enforce arc consistency
        cell = select_unassigned_cell()
        if cell is None:
            return True  # All cells are filled

        r, c = cell
        for letter in possible_values[r][c]:
            # Save the current state
            saved_grid = [row[:] for row in grid]
            saved_possible_values = [row[:] for row in possible_values]

            # Place the letter and propagate
            grid[r][c] = letter
            possible_values[r][c] = {letter}
            if backtrack():
                return True

            # Restore the state
            grid[:] = saved_grid
            possible_values[:] = saved_possible_values

        return False

    # Pre-compute feasible letters for the minor diagonal
    feasible_diag_letters = set('abcdefg')
    for r in range(7):
        feasible_diag_letters.intersection_update(possible_values[r][6-r])

    # Try each feasible letter for the minor diagonal
    for diag_letter in feasible_diag_letters:
        print(f"Trying minor diagonal letter: {diag_letter}")
        # Place diag_letter on the minor diagonal
        for r in range(7):
            grid[r][6-r] = diag_letter
            possible_values[r][6-r] = {diag_letter}

        # Try to fill the rest of the grid
        if backtrack():
            return '\n'.join(','.join(row) for row in grid)

        # Reset the grid and possible values
        grid = [list(row) for row in [
            ['g', '', 'c', 'a', 'd', '', 'e'],
            ['', 'c', '', '', '', '', ''],
            ['', '', '', '', '', 'g', ''],
            ['', '', 'b', '', 'g', 'f', ''],
            ['d', '', '', '', 'f', '', 'a'],
            ['b', '', 'g', 'f', '', 'a', ''],
            ['', '', '', '', 'a', 'd', '']
        ]]
        possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
        for r in range(7):
            for c in range(7):
                if grid[r][c] != '':
                    letter = grid[r][c]
                    possible_values[r][c] = {letter}
                    for i in range(7):
                        possible_values[r][i].discard(letter)
                        possible_values[i][c].discard(letter)

    return "No solution found"

# Initial grid
grid = [
    ['g', '', 'c', 'a', 'd', '', 'e'],
    ['', 'c', '', '', '', '', ''],
    ['', '', '', '', '', 'g', ''],
    ['', '', 'b', '', 'g', 'f', ''],
    ['d', '', '', '', 'f', '', 'a'],
    ['b', '', 'g', 'f', '', 'a', ''],
    ['', '', '', '', 'a', 'd', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(f"<<<\n{solution}\n>>>")