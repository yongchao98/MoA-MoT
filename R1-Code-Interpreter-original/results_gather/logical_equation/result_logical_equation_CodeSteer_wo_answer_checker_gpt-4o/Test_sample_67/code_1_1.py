from z3 import Int, Solver, And, Distinct

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = Int('A'), Int('B'), Int('C'), Int('D'), Int('E'), Int('F'), Int('G'), Int('H'), Int('I'), Int('J'), Int('K'), Int('L'), Int('M')

# List of possible values
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Create a solver instance
solver = Solver()

# Add constraints
solver.add(D - K == 23)
solver.add(L == 3.6 * A)
solver.add(A + H == 25)
solver.add(K == 2.5 * C)
solver.add(B > F)
solver.add(E == 2.4 * A)
solver.add(H == 1.5 * A)
solver.add(I - G == -42)
solver.add(I - E == -21)
solver.add(L == 1.5 * E)

# Ensure all variables take distinct values from the list
solver.add(Distinct(A, B, C, D, E, F, G, H, I, J, K, L, M))
solver.add(And(*[var in values for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]]))

# Check if the solution exists
if solver.check() == 'sat':
    model = solver.model()
    result = [model.evaluate(var).as_long() for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]]
    print(f"<<<{result}>>>")
else:
    print("No solution found")