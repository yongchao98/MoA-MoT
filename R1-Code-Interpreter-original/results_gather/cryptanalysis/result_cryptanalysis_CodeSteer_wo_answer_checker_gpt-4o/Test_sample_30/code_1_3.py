from constraint import Problem, AllDifferentConstraint

# Create a problem instance
problem = Problem()

# Add variables for numbers and letters
problem.addVariables("XY", range(10))  # X and Y are numbers
problem.addVariables("ZW", [chr(i) for i in range(ord('A'), ord('Z') + 1)])  # Z and W are letters

# Add constraints for numbers and letters not repeating
problem.addConstraint(AllDifferentConstraint(), "XY")
problem.addConstraint(AllDifferentConstraint(), "ZW")

# Add constraints based on feedback
# Guess 1: 38KE
problem.addConstraint(lambda X, Y: X != 3 and Y != 8, "XY")
problem.addConstraint(lambda Z, W: Z > 'K' and W > 'E', "ZW")

# Guess 2: 23BW
problem.addConstraint(lambda X, Y: X != 2 and Y != 3, "XY")
problem.addConstraint(lambda Z, W: Z != 'B' and W != 'W', "ZW")

# Guess 3: 60YQ
problem.addConstraint(lambda X, Y: (X == 6 or Y == 6) and (X != 0 and Y != 0), "XY")
problem.addConstraint(lambda Z, W: Z != 'Y' and W != 'Q', "ZW")

# Guess 4: 38LE
problem.addConstraint(lambda Z, W: Z > 'L' and W > 'E', "ZW")

# Guess 5: 18KS
problem.addConstraint(lambda X, Y: X != 1 and Y != 8, "XY")
problem.addConstraint(lambda Z, W: Z != 'K' and W != 'S', "ZW")

# Guess 6: 40NX
problem.addConstraint(lambda X, Y: (X == 4 or Y == 4) and (X != 0 and Y != 0), "XY")
problem.addConstraint(lambda Z, W: Z == 'N' and W < 'X', "ZW")

# Guess 7: 21JH
problem.addConstraint(lambda X, Y: X != 2 and Y != 1, "XY")
problem.addConstraint(lambda Z, W: Z > 'J' and W > 'H', "ZW")

# Guess 8: 73UR
problem.addConstraint(lambda X, Y: (X == 7 or Y == 7) and (X != 3 and Y != 3), "XY")
problem.addConstraint(lambda Z, W: (Z == 'U' and W != 'U') or (W == 'U' and Z != 'U'), "ZW")

# Debugging: Print intermediate solutions
solutions = problem.getSolutions()
print("Intermediate solutions:", solutions)

# Print the first valid solution
if solutions:
    solution = solutions[0]
    password = [str(solution['X']), str(solution['Y']), solution['Z'], solution['W']]
    print("Password:", password)
else:
    print("No solution found.")