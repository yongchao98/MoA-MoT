# This script calculates the position x0 based on the derived parametric solution.
# The problem's solution hinges on finding a parameter 't' from the equation:
# t^4 - 2t^3 - 3 = 0.
# One integer root of this equation is t = -1.

# The parameter value corresponding to the target point (x0, -3).
t = -1

# The position x0 is calculated from 't' using the derived relation: x0 = 3t - t^2.
# This relation is equivalent to x0 = t - 3/t^2, which was also derived in the plan.
# We use x0 = 3t - t^2 for this calculation.
x0 = 3*t - t**2

# The target vertical coordinate.
y0 = -3

print("The problem is solved by parameterizing the differential equation.")
print(f"The target condition y(x0) = {y0} corresponds to a parameter value of t = {t}.")
print("The position x0 is then found using the equation: x0 = 3*t - t^2.")
print("\nCalculation:")
print(f"x0 = 3 * ({t}) - ({t})^2")
calc_step1 = 3 * t
calc_step2 = t**2
print(f"x0 = {calc_step1} - {calc_step2}")
print(f"x0 = {x0}")
print("\n--------------------------------------------------")
print(f"The position along the trajectory where the particle reaches a vertical coordinate of {y0} is x0 = {x0}.")
print(f"The final equation with the determined position is: y({x0}) = {y0}")