from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# Define the variables and their possible values
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Add variables to the problem
problem.addVariables(letters, values)

# Add the all-different constraint
problem.addConstraint(AllDifferentConstraint())

# Add constraints based on the given equations and inequalities
problem.addConstraint(lambda D, A: D - A == 43, ('D', 'A'))
problem.addConstraint(lambda F, C: F == 1.5 * C, ('F', 'C'))
problem.addConstraint(lambda A, F: A + F == 43, ('A', 'F'))
problem.addConstraint(lambda C, L: C == 1.6 * L, ('C', 'L'))
problem.addConstraint(lambda L, J: L == 3.0 * J, ('L', 'J'))
problem.addConstraint(lambda F, J: F - J == 31, ('F', 'J'))
problem.addConstraint(lambda H, L: H > L, ('H', 'L'))
problem.addConstraint(lambda M, G: M == 4.5 * G, ('M', 'G'))
problem.addConstraint(lambda E, G: E + G == 12, ('E', 'G'))
problem.addConstraint(lambda J, K: J + K == 33, ('J', 'K'))
problem.addConstraint(lambda F, G: F == 3.6 * G, ('F', 'G'))
problem.addConstraint(lambda M, A: M > A, ('M', 'A'))

# Find a solution
solution = problem.getSolution()

# Print the solution in the required format
if solution:
    result = [solution[letter] for letter in sorted(letters)]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")