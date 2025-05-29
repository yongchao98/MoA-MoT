from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# Define the numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Add variables with their domains
problem.addVariables("ABCDEFGHIJKLM", numbers)

# Add the all-different constraint
problem.addConstraint(AllDifferentConstraint())

# Add the constraints based on the given equations and inequalities
problem.addConstraint(lambda E, H: E + H == 31, "EH")
problem.addConstraint(lambda F, K: F - K == 4, "FK")
problem.addConstraint(lambda D, K: D + K == 31, "DK")
problem.addConstraint(lambda I, L: I == 3.6 * L, "IL")
problem.addConstraint(lambda J, E: J == 1.6 * E, "JE")
problem.addConstraint(lambda C, I: C - I == -35, "CI")
problem.addConstraint(lambda A, H: A - H == -14, "AH")
problem.addConstraint(lambda I, E: I == 2.4 * E, "IE")
problem.addConstraint(lambda D, G: D + G == 33, "DG")
problem.addConstraint(lambda J, G: J > G, "JG")

# Find a solution
solution = problem.getSolution()

# Print the solution in the required format
if solution:
    result = [solution[letter] for letter in "ABCDEFGHIJKLM"]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")