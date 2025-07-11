import math

# The problem requires us to first find the value of X_0 and then compute a final expression.
# Based on the theory of boundary-value problems, we first establish a relationship between y_0 and x_0.
# The solvability condition leads to:
# y_0^4 = -x_0^6 * (7 * (10^15 - 1)) / (5 * (10^21 - 1))
# From this, we can express y_0 as a function of x_0:
# y_0(x_0) = x_0^(3/2) * [(-7 * (10^15 - 1)) / (5 * (10^21 - 1))]^(1/4)

# Next, we use the integral equation: integral from 0 to X_0 of y_0(x_0) * x_0^(p-1) dx_0 = beta
# Let K = [(-7 * (10^15 - 1)) / (5 * (10^21 - 1))]^(1/4).
# The integral becomes: K * integral from 0 to X_0 of x_0^(3/2) * x_0^5 dx_0 = beta
# K * [2/15 * X_0^(15/2)] = beta

# The value of beta is given as:
# beta = (1/1000) * (2/15) * K * 10^120
# By equating the two expressions, we can solve for X_0:
# K * (2/15) * X_0^(15/2) = (1/1000) * (2/15) * K * 10^120
# X_0^(15/2) = 10^(-3) * 10^120 = 10^117
# X_0 = (10^117)^(2/15) = 10^(234/15) = 10^15.6

# Step 1: Define and calculate X_0
exponent_X0 = 117 * (2 / 15)
X0 = 10**exponent_X0

print(f"Based on the problem's conditions, we solve for X_0.")
print(f"The equation for X_0 is: X_0^(15/2) = 10^117")
print(f"This gives X_0 = 10^({exponent_X0})")
print("-" * 30)

# Step 2: Calculate the final expression: 10^30 * X_0^2 - 10^30 * X_0 + 10
# As requested, we output each number in the final equation.
c1 = 10**30
c2 = 10**30
c3 = 10

X0_squared = X0**2
exponent_X0_squared = exponent_X0 * 2

term1 = c1 * X0_squared
term2 = c2 * X0
term3 = c3

print("The expression to calculate is: C1 * X_0^2 - C2 * X_0 + C3")
print(f"C1 = {c1:.0e}")
print(f"C2 = {c2:.0e}")
print(f"C3 = {c3}")
print(f"X_0 = 10^{exponent_X0}")
print(f"X_0^2 = 10^{exponent_X0_squared}")
print("-" * 30)

# Calculate the value of each term in the expression
exponent_term1 = 30 + exponent_X0_squared
exponent_term2 = 30 + exponent_X0

print("Calculating the terms of the final equation:")
print(f"First term: 10^30 * X_0^2 = 10^30 * 10^{exponent_X0_squared:.1f} = 10^{exponent_term1:.1f}")
print(f"Second term: 10^30 * X_0 = 10^30 * 10^{exponent_X0:.1f} = 10^{exponent_term2:.1f}")
print(f"Third term: {term3}")
print("-" * 30)

# Calculate the final result
final_result = term1 - term2 + term3

print(f"Final Result = {term1} - {term2} + {term3}")
print(f"The value of the expression is: {final_result}")