from constraint import Problem, AllDifferentConstraint

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Create a problem instance
problem = Problem()

# Add variables with their domains
problem.addVariables("ABCDEFGHIJKLM", numbers)

# Add the all-different constraint
problem.addConstraint(AllDifferentConstraint())

# Add the specific constraints
problem.addConstraint(lambda F, J: F + J == 15, ("F", "J"))
problem.addConstraint(lambda F, G: F - G == -91, ("F", "G"))
problem.addConstraint(lambda D, K: D == 3.5 * K, ("D", "K"))
problem.addConstraint(lambda E, G: E + G == 141, ("E", "G"))
problem.addConstraint(lambda A, K: A - K == 13, ("A", "K"))
problem.addConstraint(lambda F, M: F + M == 29, ("F", "M"))
problem.addConstraint(lambda D, F: D == 1.4 * F, ("D", "F"))
problem.addConstraint(lambda L, D: L == 4.0 * D, ("L", "D"))
problem.addConstraint(lambda C, J: C == 3.6 * J, ("C", "J"))
problem.addConstraint(lambda A, M: A + M == 39, ("A", "M"))
problem.addConstraint(lambda C, A: C == 2.4 * A, ("C", "A"))

# Find a solution
solution = problem.getSolution()

# Extract the solution in alphabetical order
if solution:
    result = [solution[letter] for letter in "ABCDEFGHIJKLM"]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")