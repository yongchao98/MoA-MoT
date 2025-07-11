# This script prints the derived formula for the height 'h' of the cylinder.
# The variables are:
# r: the radius of the cylinder
# theta: the inscribed angle given in the problem (in radians)
# h: the height of the cylinder

# The formula is derived from the geometry of the cylinder when its side is unrolled.
# The height 'h' forms one leg of a right triangle.
# The other leg is the arc length '2*r*theta'.
# The angle of the helical path is assumed to be 'theta'.
# tan(theta) = h / (2 * r * theta)
# Therefore, h = 2 * r * theta * tan(theta).

print("The relationship between the height 'h', the radius 'r', and the angle 'theta' is:")
print("h = 2 * r * theta * tan(theta)")
print("\nNote: For this equation to be correct, the angle 'theta' must be in radians.")
