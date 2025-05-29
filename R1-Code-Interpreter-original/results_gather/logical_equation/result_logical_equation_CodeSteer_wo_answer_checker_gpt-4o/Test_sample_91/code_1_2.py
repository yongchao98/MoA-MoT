from z3 import Int, Solver, And, Distinct

# Define variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = Int('A'), Int('B'), Int('C'), Int('D'), Int('E'), Int('F'), Int('G'), Int('H'), Int('I'), Int('J'), Int('K'), Int('L'), Int('M')

# List of possible values
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Create a solver instance
s = Solver()

# Add constraints
s.add(Distinct(A, B, C, D, E, F, G, H, I, J, K, L, M))  # All values must be distinct
s.add(And(*[v in values for v in [A, B, C, D, E, F, G, H, I, J, K, L, M]]))  # Each variable must be one of the given values

# Add the given constraints
s.add(I == 2.4 * F)
s.add(D == 4.8 * C)
s.add(D == 1.6 * F)
s.add(B == 1.6 * L)
s.add(M == 3.5 * E)
s.add(G - J == -17)
s.add(H - B == 34)
s.add(B == 3.2 * C)
s.add(H > E)
s.add(H + L == 60)
s.add(I - C == 31)
s.add(H > A)

# Check if the solution exists
if s.check() == 'sat':
    model = s.model()
    result = [model.evaluate(v).as_long() for v in [A, B, C, D, E, F, G, H, I, J, K, L, M]]
    print(f"<<<{result}>>>")
else:
    print("No solution found")