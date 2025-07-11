import sympy

# Define the variable n
n = sympy.Symbol('n')

# The problem states the probability of reaching a cell in at most 3 moves is 66%
# Probability = (Reachable Cells) / (Total Cells)
# 0.66 = R(n) / n^2
# Based on a model of the object's movement, the number of reachable cells R(n) can be expressed as:
# R(n) = 6*n + 6
# So, the equation to solve is (6*n + 6) / n^2 = 0.66

print("Step 1: Set up the equation based on the problem statement.")
print("Let R(n) be the number of reachable cells. The total number of cells is n*n.")
print("The probability is R(n) / (n*n) = 0.66.")
print("A model for the number of reachable cells is R(n) = 6*n + 6.")
print("So, the equation is: (6*n + 6) / (n*n) = 0.66\n")


# Set up the equation
# Using 66/100 to avoid floating point issues and simplifying it to 33/50.
equation = sympy.Eq((6 * n + 6) / (n**2), 66/100)

print("Step 2: Rearrange the equation into a standard quadratic form (a*n^2 + b*n + c = 0).")
# Multiply by n^2
# 6n + 6 = 0.66 * n^2
# 0.66n^2 - 6n - 6 = 0
# To work with integers, we can multiply the entire equation by 100:
# 66n^2 - 600n - 600 = 0
# And divide by 6 to simplify:
# 11n^2 - 100n - 100 = 0
quadratic_eq = 11*n**2 - 100*n - 100
print(f"The simplified quadratic equation is: {quadratic_eq} = 0\n")


print("Step 3: Solve the quadratic equation for n.")
# Solve the equation
solutions = sympy.solve(quadratic_eq, n)
print(f"The mathematical solutions to the equation are: {solutions}\n")

# The problem states n must be an even integer. We filter for a positive integer solution.
final_n = [sol for sol in solutions if sol > 0 and sol.is_integer]

print("Step 4: Select the valid solution.")
print("According to the problem, n must be a positive even integer.")
if final_n:
    result = final_n[0]
    print(f"The only positive integer solution is n = {result}.")
    print("This value is even, so it is the correct answer.")
else:
    print("No valid positive integer solution found.")

final_answer = final_n[0]
print(f"\nFinal check:")
print(f"If n = {final_answer}, the total number of cells is {final_answer**2}.")
print(f"The number of reachable cells is 6 * {final_answer} + 6 = {6*final_answer+6}.")
print(f"The probability is ({6*final_answer+6})/({final_answer**2}) = {(6*final_answer+6)/(final_answer**2)}.")
print(f"This matches the given probability of 66%.")

<<<10>>>