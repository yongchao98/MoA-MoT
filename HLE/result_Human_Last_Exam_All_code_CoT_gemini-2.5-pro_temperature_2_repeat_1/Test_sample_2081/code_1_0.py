import math

# The problem is to find the radius R of the sphere of initial conditions
# (x0, y0, z0) for which the nonlinear system has a solution.
# The analysis leads to the equation for the radius R, which is
# R^2 = 0.5 * e^T * (e^T + 1).

# Given T = ln(10^34), so e^T = 10^34.
# We will now calculate R based on this formula.

# The final equation is R = sqrt(term1 * term2 * term3)
term1 = 0.5
term2 = 10**34
term3 = 10**34 + 1

# Note that due to standard floating-point precision, term3 will be numerically
# indistinguishable from term2. However, the symbolic formula is correct.
# We calculate R based on these values.

R_squared = term1 * term2 * term3
R = math.sqrt(R_squared)

print("The final equation for R is of the form: R = sqrt(c1 * c2 * c3)")
print(f"The value for c1 is: {term1}")
print(f"The value for c2 is: {term2}")
print(f"The value for c3 is: {term3}")
print(f"The calculated value for R-squared is: {R_squared}")
print(f"The final calculated value for R is: {R}")