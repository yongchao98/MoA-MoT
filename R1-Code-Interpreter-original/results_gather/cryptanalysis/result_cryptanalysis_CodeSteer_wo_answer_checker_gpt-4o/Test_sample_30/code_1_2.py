from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# Add variables for numbers and letters
problem.addVariables("XY", range(10))  # X and Y are numbers
problem.addVariables("ZW", [chr(i) for i in range(ord('A'), ord('Z') + 1)])  # Z and W are letters

# Add constraints based on the feedback
problem.addConstraint(AllDifferentConstraint(), "XY")
problem.addConstraint(AllDifferentConstraint(), "ZW")

# Feedback constraints
problem.addConstraint(lambda X, Y: X != 3 and Y != 8, "XY")
problem.addConstraint(lambda Z, W: Z > 'K' and W > 'E', "ZW")
problem.addConstraint(lambda X, Y: X != 2 and Y != 3, "XY")
problem.addConstraint(lambda Z, W: Z != 'B' and W != 'W', "ZW")
problem.addConstraint(lambda X, Y: (X == 6 or Y == 6) and (X != 0 and Y != 0), "XY")
problem.addConstraint(lambda Z, W: Z != 'Y' and W != 'Q', "ZW")
problem.addConstraint(lambda X, Y: X != 1 and Y != 8, "XY")
problem.addConstraint(lambda Z, W: Z != 'K' and W != 'S', "ZW")
problem.addConstraint(lambda X, Y: (X == 4 or Y == 4) and (X != 0 and Y != 0), "XY")
problem.addConstraint(lambda Z, W: Z == 'N' and W < 'X', "ZW")
problem.addConstraint(lambda X, Y: X != 2 and Y != 1, "XY")
problem.addConstraint(lambda Z, W: Z > 'J' and W > 'H', "ZW")
problem.addConstraint(lambda X, Y: (X == 7 or Y == 7) and (X != 3 and Y != 3), "XY")
problem.addConstraint(lambda Z, W: (Z == 'U' or W == 'U') and (Z != 'R' and W != 'R'), "ZW")

# Adjust constraints for correct numbers/letters in wrong positions
problem.addConstraint(lambda X, Y: (X == 6 and Y != 6) or (Y == 6 and X != 6), "XY")
problem.addConstraint(lambda X, Y: (X == 4 and Y != 4) or (Y == 4 and X != 4), "XY")
problem.addConstraint(lambda X, Y: (X == 7 and Y != 7) or (Y == 7 and X != 7), "XY")
problem.addConstraint(lambda Z, W: (Z == 'U' and W != 'U') or (W == 'U' and Z != 'U'), "ZW")

# Find solutions
solutions = problem.getSolutions()

# Print the first valid solution
if solutions:
    solution = solutions[0]
    password = [str(solution['X']), str(solution['Y']), solution['Z'], solution['W']]
    print(password)
else:
    print("No solution found.")