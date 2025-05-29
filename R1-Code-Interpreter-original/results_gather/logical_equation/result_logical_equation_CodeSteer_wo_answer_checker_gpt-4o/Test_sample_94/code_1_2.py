from z3 import Int, Solver, And, Distinct

# Define variables for each letter
A = Int('A')
B = Int('B')
C = Int('C')
D = Int('D')
E = Int('E')
F = Int('F')
G = Int('G')
H = Int('H')
I = Int('I')
J = Int('J')
K = Int('K')
L = Int('L')
M = Int('M')

# List of all variables
variables = [A, B, C, D, E, F, G, H, I, J, K, L, M]

# Possible values for the variables
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Create a solver instance
solver = Solver()

# Add constraints
solver.add(Distinct(variables))  # All variables must have distinct values
solver.add(And([var >= 2 for var in variables]))  # All values must be at least 2
solver.add(And([var <= 96 for var in variables]))  # All values must be at most 96

# Add specific constraints
solver.add(G + I == 55)
solver.add(E > D)
solver.add(G + L == 57)
solver.add(G - E == 5)
solver.add(C + H == 27)
solver.add(A - I == 5)
solver.add(F - I == 91)
solver.add(B - C == -9)
solver.add(K > C)
solver.add(H + L == 10)
solver.add(J - D == 14)
solver.add(B - F == -81)

# Add constraints to ensure each variable takes one of the given values
solver.add(And([Or([var == val for val in values]) for var in variables]))

# Check if the problem is solvable
if solver.check() == 'sat':
    model = solver.model()
    result = [model.evaluate(var).as_long() for var in variables]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")