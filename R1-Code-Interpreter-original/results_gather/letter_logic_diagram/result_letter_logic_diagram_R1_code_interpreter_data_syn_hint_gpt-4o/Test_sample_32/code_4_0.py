from ortools.sat.python import cp_model

# Create the model
model = cp_model.CpModel()

# Create the variables
grid = {}
for i in range(7):
    for j in range(7):
        grid[(i, j)] = model.NewIntVar(0, 6, f'grid_{i}_{j}')

# Pre-filled grid
prefilled = [
    [None, 1, None, 2, 3, 5, 0],
    [None, None, 2, None, None, None, 6],
    [None, 2, None, None, None, 6, None],
    [2, 3, None, 0, None, None, 4],
    [3, None, 0, None, 1, 4, None],
    [5, 0, None, None, None, None, None],
    [0, None, 1, 4, 2, None, 5]
]

# Add constraints for pre-filled cells
for i in range(7):
    for j in range(7):
        if prefilled[i][j] is not None:
            model.Add(grid[(i, j)] == prefilled[i][j])

# Add constraints for rows and columns
for i in range(7):
    model.AddAllDifferent([grid[(i, j)] for j in range(7)])
    model.AddAllDifferent([grid[(j, i)] for j in range(7)])

# Add constraint for the minor diagonal
minor_diagonal_value = model.NewIntVar(0, 6, 'minor_diagonal_value')
for i in range(7):
    model.Add(grid[(i, 6 - i)] == minor_diagonal_value)

# Solve the model
solver = cp_model.CpSolver()
solver.Solve(model)

# Print the solution
solution = []
for i in range(7):
    row = []
    for j in range(7):
        row.append(chr(solver.Value(grid[(i, j)]) + ord('a')))
    solution.append(",".join(row))

for row in solution:
    print(row)