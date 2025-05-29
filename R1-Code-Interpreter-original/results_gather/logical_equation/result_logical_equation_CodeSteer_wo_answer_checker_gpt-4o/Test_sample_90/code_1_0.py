from constraint import Problem

def solve_with_constraint_library():
    # Create a problem instance
    problem = Problem()

    # Define the variables and their domains
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    variables = "ABCDEFGHIJKLM"
    for var in variables:
        problem.addVariable(var, numbers)

    # Add constraints
    problem.addConstraint(lambda E, K: E + K == 81, ("E", "K"))
    problem.addConstraint(lambda K, J: K > J, ("K", "J"))
    problem.addConstraint(lambda E, G: E == 1.5 * G, ("E", "G"))
    problem.addConstraint(lambda B, J: B + J == 5, ("B", "J"))
    problem.addConstraint(lambda M, F: M - F == -6, ("M", "F"))
    problem.addConstraint(lambda K, M: K == 4.5 * M, ("K", "M"))
    problem.addConstraint(lambda M, B: M > B, ("M", "B"))
    problem.addConstraint(lambda F, C: F == 3.2 * C, ("F", "C"))
    problem.addConstraint(lambda C, E: C - E == -31, ("C", "E"))
    problem.addConstraint(lambda L, C: L == 1.4 * C, ("L", "C"))
    problem.addConstraint(lambda A, M: A + M == 38, ("A", "M"))
    problem.addConstraint(lambda F, I: F + I == 112, ("F", "I"))

    # Solve the problem
    solution = problem.getSolution()

    if solution:
        # Extract the solution in alphabetical order
        result = [solution[var] for var in sorted(variables)]
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

solve_with_constraint_library()