from z3 import Int, Solver, And, Distinct

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Create a solver instance
solver = Solver()

# Define variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = Int('A'), Int('B'), Int('C'), Int('D'), Int('E'), Int('F'), Int('G'), Int('H'), Int('I'), Int('J'), Int('K'), Int('L'), Int('M')

# Add constraints for the domain of each variable
solver.add(Distinct(A, B, C, D, E, F, G, H, I, J, K, L, M))
for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]:
    solver.add(Or([var == n for n in numbers]))

# Add the given constraints
solver.add(C - E == -21)
solver.add(J == 4 * C)
solver.add(H + K == 23)
solver.add(D + L == 17)
solver.add(C - A == -26)
solver.add(G + M == 8)
solver.add(F + G == 13)
solver.add(H - G == 4)
solver.add(A - C == 26)
solver.add(I - J == -68)
solver.add(F + M == 15)

# Solve the problem
if solver.check() == 'sat':
    model = solver.model()
    solution = [model.evaluate(var).as_long() for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]]
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")