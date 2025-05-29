from z3 import Int, Solver, And, Distinct

# Define the variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = Int('A'), Int('B'), Int('C'), Int('D'), Int('E'), Int('F'), Int('G'), Int('H'), Int('I'), Int('J'), Int('K'), Int('L'), Int('M')

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Create a solver instance
solver = Solver()

# Add constraints for the equations and inequalities
solver.add(D - A == 43)
solver.add(F == 1.5 * C)
solver.add(A + F == 43)
solver.add(C == 1.6 * L)
solver.add(L == 3.0 * J)
solver.add(F - J == 31)
solver.add(H > L)
solver.add(M == 4.5 * G)
solver.add(E + G == 12)
solver.add(J + K == 33)
solver.add(F == 3.6 * G)
solver.add(M > A)

# Add constraints for the possible values
solver.add(Distinct(A, B, C, D, E, F, G, H, I, J, K, L, M))
for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]:
    solver.add(Or([var == v for v in values]))

# Check if a solution exists
if solver.check() == 'sat':
    model = solver.model()
    result = [model.evaluate(var).as_long() for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")