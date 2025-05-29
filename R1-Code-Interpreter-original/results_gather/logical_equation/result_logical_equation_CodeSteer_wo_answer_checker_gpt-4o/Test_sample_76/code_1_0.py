from constraint import Problem

# Create a problem instance
problem = Problem()

# Define the variables and their possible values
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Add variables to the problem
problem.addVariables(letters, numbers)

# Add constraints based on the given equations
problem.addConstraint(lambda I, E: I - E == -4, ('I', 'E'))
problem.addConstraint(lambda F, K: F == 2.5 * K, ('F', 'K'))
problem.addConstraint(lambda I, B: I == 1.5 * B, ('I', 'B'))
problem.addConstraint(lambda L, A: L - A == 33, ('L', 'A'))
problem.addConstraint(lambda E, H: E == 2.8 * H, ('E', 'H'))
problem.addConstraint(lambda C, M: C + M == 111, ('C', 'M'))
problem.addConstraint(lambda L, C: L - C == -60, ('L', 'C'))
problem.addConstraint(lambda H, F: H == 2.0 * F, ('H', 'F'))
problem.addConstraint(lambda F, L: F + L == 41, ('F', 'L'))
problem.addConstraint(lambda K, E: K - E == -26, ('K', 'E'))
problem.addConstraint(lambda A, F: A + F == 8, ('A', 'F'))

# Ensure all variables have unique values
problem.addConstraint(lambda *args: len(set(args)) == len(args), letters)

# Find a solution
solution = problem.getSolution()

# Print the solution in the required format
if solution:
    result = [solution[letter] for letter in sorted(letters)]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")