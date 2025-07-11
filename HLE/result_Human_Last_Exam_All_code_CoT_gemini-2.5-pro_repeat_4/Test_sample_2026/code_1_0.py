import math

# From the analysis of the problem, the possible initial values (x0, y0, z0)
# lie on the surface of a sphere defined by the equation x0^2 + y0^2 + z0^2 = R^2.
# The radius squared R^2 is derived from the problem's conditions.
r_squared = 10**25 / 2

# The area requested is the surface area of this sphere.
# The formula for the surface area of a sphere is A = 4 * pi * R^2.
# We will now calculate this area.

# The final equation for the area is Area = 4 * pi * r_squared.
# As requested, we will output each number in this equation.
four = 4
pi_val = math.pi

print("The set of initial values (x0, y0, z0) forms a sphere.")
print(f"The equation of the sphere is: x0^2 + y0^2 + z0^2 = {r_squared}")
print("\nThe area of this surface is calculated using the formula: Area = 4 * pi * R^2")
print("\nThe numbers in this final equation are:")
print(f"The number 4: {four}")
print(f"The value of pi: {pi_val}")
print(f"The value of R^2: {r_squared}")

# Calculate the final surface area
surface_area = four * pi_val * r_squared

print(f"\nThe calculated surface area is: {surface_area}")
