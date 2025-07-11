import math

# Define the parameters from the problem statement
A = 10**10
B = 1e-5 - 1
T = 10**20

# As derived from the differential equation and its boundary condition,
# the solution X_0(t) is a constant given by the formula:
# X_0(t) = A / (B + 1).
# The question asks for the value of X_0(t) at t=T.

# Define the components of the final equation
numerator = A
denominator = B + 1

# Calculate the result
result = numerator / denominator

# Output the equation with the given numbers and the final result,
# as requested.
print("The value of X_0(t) is constant and calculated using the formula: X_0(T) = A / (B + 1)")
print(f"where T = {T:.0e}, A = {numerator:.0e}, and B + 1 = {denominator:.0e}")
print("Substituting the numbers into the equation:")
print(f"X_0({T:.0e}) = {numerator:.0e} / {denominator:.0e}")
print(f"Result: X_0({T:.0e}) = {result:.0e}")