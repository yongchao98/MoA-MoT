from ortools.sat.python import cp_model

def solve_puzzle():
    model = cp_model.CpModel()
    n = 7
    letters = 'abcdefg'
    letter_to_num = {letter: i for i, letter in enumerate(letters)}
    num_to_letter = {i: letter for i, letter in enumerate(letters)}

    # Create variables
    grid = [[model.NewIntVar(0, 6, f'cell_{i}_{j}') for j in range(n)] for i in range(n)]

    # Pre-filled grid
    pre_filled = [
        ['d', 'a', '', '', 'g', '', ''],
        ['', '', '', 'g', 'e', '', ''],
        ['f', '', 'g', '', 'c', '', ''],
        ['', 'g', '', '', '', 'a', ''],
        ['g', 'e', '', '', '', 'f', ''],
        ['e', '', '', '', '', 'b', ''],
        ['', 'd', '', 'f', 'b', '', 'e']
    ]

    # Apply pre-filled constraints
    for i in range(n):
        for j in range(n):
            if pre_filled[i][j]:
                model.Add(grid[i][j] == letter_to_num[pre_filled[i][j]])

    # Row and column constraints
    for i in range(n):
        model.AddAllDifferent([grid[i][j] for j in range(n)])
        model.AddAllDifferent([grid[j][i] for j in range(n)])

    # Minor diagonal constraint
    minor_diagonal_indices = [(i, n - 1 - i) for i in range(n)]
    minor_diagonal_letter = model.NewIntVar(0, 6, 'minor_diagonal_letter')
    for i, j in minor_diagonal_indices:
        model.Add(grid[i][j] == minor_diagonal_letter)

    # Solve the model
    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    if status == cp_model.OPTIMAL:
        solution = []
        for i in range(n):
            row = [num_to_letter[solver.Value(grid[i][j])] for j in range(n)]
            solution.append(','.join(row))
        return solution
    else:
        return None

# Get the solution
solution = solve_puzzle()

# Print the solution
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")