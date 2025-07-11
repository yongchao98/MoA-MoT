x1 = 6
x2 = 9
x3 = 2
x4 = 8
x5 = 7
x6 = 4
x7 = 0
x8 = 6
x9 = 9
x10 = 0
x11 = 4
x12 = 2

# The equation is: X1*X2 + (10*X3+X4 - X5*X6) + X7*X4 = X8*X9
# Let's verify it
lhs = x1 * x2 + (10 * x3 + x4 - x5 * x6) + x7 * x4
rhs = x8 * x9

# The problem requires printing the numbers in the final equation.
# Let's construct the string for the equation
equation_str = f"{x1} * {x2} + ({x3} * 10 + {x4} - {x5} * {x6}) + {x7} * {x4} = {x8} * {x9}"

print("The solved equation is:")
print(equation_str)
print(f"{lhs} = {rhs}")
print("\nThe value of X12 is:")
print(x12)