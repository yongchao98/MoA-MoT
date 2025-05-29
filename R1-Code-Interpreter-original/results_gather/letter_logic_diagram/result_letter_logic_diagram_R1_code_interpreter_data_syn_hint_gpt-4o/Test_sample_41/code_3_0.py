from ortools.sat.python import cp_model

# Create the model
model = cp_model.CpModel()

# Create the variables
grid = {}
letters = 'abcdefg'
for r in range(7):
    for c in range(7):
        if (r, c) in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
            grid[(r, c)] = model.NewIntVar(0, 6, f'grid_{r}_{c}')
        else:
            grid[(r, c)] = model.NewIntVar(0, 6, f'grid_{r}_{c}')

# Add the constraints for rows and columns
for i in range(7):
    model.AddAllDifferent([grid[(i, j)] for j in range(7)])
    model.AddAllDifferent([grid[(j, i)] for j in range(7)])

# Add the constraint for the minor diagonal
minor_diagonal_value = model.NewIntVar(0, 6, 'minor_diagonal_value')
for r, c in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
    model.Add(grid[(r, c)] == minor_diagonal_value)

# Pre-fill the known values
pre_filled = [
    (0, 1, 'a'), (0, 2, 'g'),
    (1, 0, 'a'), (1, 1, 'g'), (1, 3, 'f'), (1, 4, 'b'), (1, 6, 'd'),
    (2, 1, 'c'),
    (3, 1, 'f'), (3, 2, 'b'), (3, 4, 'd'),
    (4, 5, 'g'),
    (5, 1, 'e'), (5, 4, 'g'), (5, 5, 'c'),
    (6, 2, 'a'), (6, 3, 'g'), (6, 4, 'c'), (6, 6, 'b')
]

for r, c, letter in pre_filled:
    model.Add(grid[(r, c)] == letters.index(letter))

# Create the solver and solve
solver = cp_model.CpSolver()
status = solver.Solve(model)

# Output the solution
if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
    solution = []
    for r in range(7):
        row = []
        for c in range(7):
            row.append(letters[solver.Value(grid[(r, c)])])
        solution.append(",".join(row))
    for row in solution:
        print(row)