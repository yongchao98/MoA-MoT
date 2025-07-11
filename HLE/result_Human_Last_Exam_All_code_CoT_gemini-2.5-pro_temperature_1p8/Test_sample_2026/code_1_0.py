import math

# As derived in the explanation, the condition on the initial values (x0, y0, z0)
# is that they must lie on a sphere defined by the equation:
# x0^2 + y0^2 + z0^2 = R^2
# where the radius squared is:
R_squared = 10**25 / 2.0

# The problem asks for the area of this surface.
# The area of a sphere is calculated using the formula A = 4 * pi * R^2.
four = 4.0
pi_val = math.pi
area = four * pi_val * R_squared

# The following print statements show the calculation step-by-step
# and output each number used in the final equation for the area.
print("The values (x0, y0, z0) must lie on a sphere with a radius squared R^2, calculated as:")
print(f"R^2 = 10^25 / 2 = {R_squared:.1e}")
print("\nThe surface area of this sphere is calculated with the formula: Area = 4 * pi * R^2.")
print("The final calculation, showing each number in the equation, is:")
print(f"Area = {four} * {pi_val} * {R_squared:.1e}")
print(f"Area = {area}")
print(f"\nThe final area is approximately {area:.4e}.")
