from fractions import Fraction

# This script calculates x(1) for the initial value problem
# x'(t) = (t-1) * x^2(t), with x(0) = -8.
# The analytical solution leads to the equation for x(1):
# -1/x(1) = (1)^2/2 - 1 + 1/8

print("Solving for x(1) from the equation: -1/x(1) = (1)^2/2 - 1 + 1/8")

# Define the numbers in the equation using the Fraction class for exact arithmetic
val1 = Fraction(1**2, 2)
val2 = Fraction(-1, 1)
val3 = Fraction(1, 8)

# Calculate the right-hand side (RHS) of the equation
rhs = val1 + val2 + val3

# Calculate the final result for x(1), which is -1 / RHS
x1 = -1 / rhs

print(f"The numbers on the right-hand side of the equation are: {val1}, {val2}, and {val3}")
print(f"The sum on the right-hand side is: {rhs}")
print(f"The equation becomes: -1/x(1) = {rhs}")
print(f"Therefore, x(1) = -1 / ({rhs})")
print(f"The final exact result is: x(1) = {x1}")
print(f"As a decimal, the result is approximately: {float(x1)}")