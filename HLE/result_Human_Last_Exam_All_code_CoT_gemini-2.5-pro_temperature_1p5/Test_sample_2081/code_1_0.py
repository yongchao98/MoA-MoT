import math

# From the problem statement, T = ln(10^34), which means e^T = 10^34.
exp_T = 10**34

# The set of initial values (x₀, y₀, z₀) for which solutions exist
# forms a sphere described by the equation:
# 1*x₀² + 1*y₀² + 1*z₀² = R²
# The radius squared R² is calculated using the derived formula:
# R² = 0.5 * e^T * (e^T + 1)

# Calculate R² using the simplified and numerically stable formula.
R_squared = 0.5 * exp_T * (exp_T + 1)

# Calculate the radius R.
R = math.sqrt(R_squared)

# The "final equation" is the equation of the sphere. We output its coefficients and right-hand side.
# The coefficients of x₀², y₀², and z₀² are all 1.
print(f"The equation for the surface of initial values is: 1*x₀² + 1*y₀² + 1*z₀² = {R_squared}")

# Finally, we output the value of R as requested by the problem.
print(f"The radius R is: {R}")