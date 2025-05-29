from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# Define the variables and their domains
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
problem.addVariables(letters, numbers)

# Add the constraint that all variables must have different values
problem.addConstraint(AllDifferentConstraint())

# Add the constraints based on the given conditions
problem.addConstraint(lambda B, I: B + I == 156, ('B', 'I'))
problem.addConstraint(lambda B, D: B + D == 123, ('B', 'D'))
problem.addConstraint(lambda D, A: D - A == -29, ('D', 'A'))
problem.addConstraint(lambda G, D: G == 2 * D, ('G', 'D'))
problem.addConstraint(lambda B, C: B == 1.5 * C, ('B', 'C'))
problem.addConstraint(lambda A, D: A > D, ('A', 'D'))
problem.addConstraint(lambda A, I: A - I == -4, ('A', 'I'))
problem.addConstraint(lambda C, F: C - F == 71, ('C', 'F'))

# Find a solution
solution = problem.getSolution()

# Extract the solution in alphabetical order
if solution:
    result = [solution[letter] for letter in sorted(letters)]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")