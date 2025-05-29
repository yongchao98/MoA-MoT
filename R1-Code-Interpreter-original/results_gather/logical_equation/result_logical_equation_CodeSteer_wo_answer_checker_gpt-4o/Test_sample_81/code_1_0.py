from constraint import Problem, AllDifferentConstraint

# Define the problem
problem = Problem()

# Define the variables and their possible values
variables = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Add variables to the problem
for var in variables:
    problem.addVariable(var, possible_values)

# Add the all-different constraint
problem.addConstraint(AllDifferentConstraint())

# Add the specific constraints
problem.addConstraint(lambda F, G: F + G == 46, ('F', 'G'))
problem.addConstraint(lambda J, M: J - M == 80, ('J', 'M'))
problem.addConstraint(lambda J, F: J - F == 86, ('J', 'F'))
problem.addConstraint(lambda H, L: H - L == 26, ('H', 'L'))
problem.addConstraint(lambda B, M: B + M == 31, ('B', 'M'))
problem.addConstraint(lambda B, F: B == 1.5 * F, ('B', 'F'))
problem.addConstraint(lambda C, A: C == 4.0 * A, ('C', 'A'))
problem.addConstraint(lambda L, M: L == 1.5 * M, ('L', 'M'))
problem.addConstraint(lambda G, F: G == 3.6 * F, ('G', 'F'))
problem.addConstraint(lambda D, E: D == 1.5 * E, ('D', 'E'))
problem.addConstraint(lambda L, I: L - I == -21, ('L', 'I'))

# Find a solution
solution = problem.getSolution()

# Extract the values in alphabetical order and output the result
if solution:
    result = [solution[var] for var in sorted(variables)]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")