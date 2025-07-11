import math

# The problem is to find the shape of a given volume that maximizes gravitational attraction at a point.
# This classic problem results in a shape defined in spherical coordinates by r = R_max * sqrt(cos(theta)).
# The point of maximum attraction, A, is at the origin.
# The volume (V) of this shape is related to its maximum radius (R_max) by the formula:
# V = (4 * pi / 15) * R_max^3

# We are given that the volume is 1 cubic meter.
V = 1.0

# We need to solve for R_max, which represents the furthest point on the surface from point A.
# R_max^3 = V * 15 / (4 * pi)
# R_max = (V * 15 / (4 * pi))^(1/3)

# Calculate R_max
r_max = (V * 15 / (4 * math.pi))**(1/3)

print("The formula for the furthest point R_max is (15 / (4 * pi))^(1/3).")
print(f"Given V = 1 m^3:")
print(f"R_max = (15 / (4 * {math.pi}))^(1/3)")
print(f"R_max = ({15 / (4 * math.pi)})^(1/3)")
print(f"The distance to the furthest point on the surface is: {r_max} meters.")
<<<1.0607953281144047>>>