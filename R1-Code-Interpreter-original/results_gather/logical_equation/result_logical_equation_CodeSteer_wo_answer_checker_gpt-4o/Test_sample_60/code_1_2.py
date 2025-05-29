from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# List of numbers to assign
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Add variables with their domains
problem.addVariables("ABCDEFGHIJKLM", numbers)

# Add the all-different constraint
problem.addConstraint(AllDifferentConstraint())

# Add the specific constraints
problem.addConstraint(lambda A, G: A - G == -13, "AG")
problem.addConstraint(lambda B, C: B + C == 15, "BC")
problem.addConstraint(lambda H, B: H == 4.8 * B, "HB")
problem.addConstraint(lambda I, B: I - B == -3, "IB")
problem.addConstraint(lambda L, J: L - J == -5, "LJ")
problem.addConstraint(lambda B, A: B > A, "BA")
problem.addConstraint(lambda I, F: I == 2.0 * F, "IF")
problem.addConstraint(lambda C, H: C - H == -14, "CH")
problem.addConstraint(lambda B, I: B - I == 3, "BI")
problem.addConstraint(lambda A, F: A - F == 2, "AF")

# Solve the problem
solution = problem.getSolutions()

# If a solution is found, print it
if solution:
    # Extract the first solution and sort it by variable name
    first_solution = solution[0]
    result = [first_solution[var] for var in "ABCDEFGHIJKLM"]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")