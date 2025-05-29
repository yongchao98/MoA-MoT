from z3 import Int, Solver, Distinct

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = Int('A'), Int('B'), Int('C'), Int('D'), Int('E'), Int('F'), Int('G'), Int('H'), Int('I'), Int('J'), Int('K')

# List of possible values
values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Create a solver instance
s = Solver()

# Add constraints for each variable to be one of the given numbers
s.add(A in values, B in values, C in values, D in values, E in values, F in values, G in values, H in values, I in values, J in values, K in values)

# Add the distinct constraint
s.add(Distinct(A, B, C, D, E, F, G, H, I, J, K))

# Add the given equations and inequalities
s.add(I == 4 * C)
s.add(F + I == 83)
s.add(K - E == 82)
s.add(E + F == 42)
s.add(A + H == 225)
s.add(B + K == 130)
s.add(A > B)
s.add(J == 3 * H)
s.add(B == 3 * F)
s.add(C + F == 23)
s.add(A > K)

# Check if the constraints are satisfiable
if s.check() == 'sat':
    # Get the model (solution)
    model = s.model()
    # Extract the values in alphabetical order
    solution = [model[A], model[B], model[C], model[D], model[E], model[F], model[G], model[H], model[I], model[J], model[K]]
    # Print the solution
    print(f"<<<{solution}>>>")
else:
    print("No solution found")