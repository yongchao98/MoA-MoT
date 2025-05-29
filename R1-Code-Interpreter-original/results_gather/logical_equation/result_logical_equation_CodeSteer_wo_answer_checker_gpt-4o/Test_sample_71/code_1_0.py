from constraint import Problem, AllDifferentConstraint

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Create a problem instance
problem = Problem()

# Add variables
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
problem.addVariables(letters, numbers)

# Add the all-different constraint
problem.addConstraint(AllDifferentConstraint())

# Add constraints based on the equations
problem.addConstraint(lambda B, F: B - F == -21, ('B', 'F'))
problem.addConstraint(lambda F, L: F + L == 73, ('F', 'L'))
problem.addConstraint(lambda F, A: F - A == -22, ('F', 'A'))
problem.addConstraint(lambda K, L: K - L == -42, ('K', 'L'))
problem.addConstraint(lambda B, C: B + C == 43, ('B', 'C'))
problem.addConstraint(lambda H, L: H + L == 60, ('H', 'L'))
problem.addConstraint(lambda H, M: H + M == 17, ('H', 'M'))
problem.addConstraint(lambda C, H: C == 2.4 * H, ('C', 'H'))
problem.addConstraint(lambda D, K: D + K == 4, ('D', 'K'))
problem.addConstraint(lambda G, M: G == 2.5 * M, ('G', 'M'))

# Find a solution
solution = problem.getSolution()

# Print the solution in the required format
if solution:
    result = [solution[letter] for letter in sorted(letters)]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")