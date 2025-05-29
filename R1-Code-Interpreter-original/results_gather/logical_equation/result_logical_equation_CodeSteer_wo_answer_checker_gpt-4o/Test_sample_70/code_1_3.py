from constraint import Problem

# Create a problem instance
problem = Problem()

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Add variables with their domains
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
for letter in letters:
    problem.addVariable(letter, numbers)

# Add constraints
problem.addConstraint(lambda M, L: M == 2.4 * L, ('M', 'L'))
problem.addConstraint(lambda F, L: F == 3.0 * L, ('F', 'L'))
problem.addConstraint(lambda E, I: E + I == 8, ('E', 'I'))
problem.addConstraint(lambda K, H: K - H == -8, ('K', 'H'))
problem.addConstraint(lambda H, I: H > I, ('H', 'I'))
problem.addConstraint(lambda L, M: L + M == 51, ('L', 'M'))
problem.addConstraint(lambda C, J: C + J == 9, ('C', 'J'))
problem.addConstraint(lambda G, I: G == 2.0 * I, ('G', 'I'))
problem.addConstraint(lambda A, D: A - D == 49, ('A', 'D'))
problem.addConstraint(lambda H, K: H == 1.5 * K, ('H', 'K'))

# Solve the problem
solution = problem.getSolution()

# Print the solution in the required format
if solution:
    result = [solution[letter] for letter in letters]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")