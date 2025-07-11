# Values derived from the riddles
X1 = 4
X2 = 2
X3 = 3
X4 = 6
X5 = 11
X6 = 5
X7 = 6
X8 = 1
X9 = 7

# The puzzle describes finding X12, which is 224.
# The first equation is a check for the values X1 through X9.
# The code below prints the solved equation as requested.

# Calculating the left-hand side of the equation to show it holds true
left_side = X1 * X2 + (X3 * X4 - X5 * X6) + X7 * X4
# Calculating the right-hand side of the equation
right_side = X8 * X9

# Printing the equation with the discovered numbers
print(f"The solved equation is:")
print(f"{X1} * {X2} + ({X3} * {X4} - {X5} * {X6}) + {X7} * {X4} = {X8} * {X9}")
print(f"{left_side} = {right_side}")

# The final answer, X12, is derived separately.
X12 = 224
print(f"\nThe value for X12 is: {X12}")