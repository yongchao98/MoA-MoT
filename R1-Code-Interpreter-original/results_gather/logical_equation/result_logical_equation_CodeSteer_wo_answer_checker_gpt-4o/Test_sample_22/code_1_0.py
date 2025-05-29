from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# Possible values for the letters
values = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Add variables with their possible values
problem.addVariables(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'], values)

# Add the constraint that all variables must have different values
problem.addConstraint(AllDifferentConstraint())

# Add the constraints based on the problem description
problem.addConstraint(lambda B, G: B == 2 * G, ['B', 'G'])
problem.addConstraint(lambda I, F: I == 2.5 * F, ['I', 'F'])
problem.addConstraint(lambda B, I: B + I == 86, ['B', 'I'])
problem.addConstraint(lambda A, G: A + G == 12, ['A', 'G'])
problem.addConstraint(lambda I, F: I > F, ['I', 'F'])
problem.addConstraint(lambda F, I: F + I == 112, ['F', 'I'])
problem.addConstraint(lambda C, H: C + H == 140, ['C', 'H'])
problem.addConstraint(lambda A, H: A - H == -11, ['A', 'H'])

# Find a solution
solution = problem.getSolution()

# Extract the solution in alphabetical order
result = [solution[letter] for letter in sorted(solution.keys())]

# Print the solution in the required format
print(f"<<<{result}>>>")