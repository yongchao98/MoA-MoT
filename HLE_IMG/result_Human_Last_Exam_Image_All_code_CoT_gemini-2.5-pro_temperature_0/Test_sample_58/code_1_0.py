import math

# The problem asks for the height 'h' in terms of the radius 'r' and angle 'theta'.
# Based on geometric analysis, we can establish a relationship.

# We form a right-angled triangle in the vertical plane with sides:
# - The diameter of the cylinder, which is 2*r.
# - The height of the cylinder, h.

# The tangent of the angle (let's call it alpha) in this triangle is tan(alpha) = h / (2*r).
# The problem implies a geometric analogy where this angle alpha is equal to the given angle theta.
# So, tan(theta) = h / (2*r).

# Solving for h, we get: h = 2 * r * tan(theta).

# The following code prints this final equation.
# The number '2' is explicitly part of the equation as requested.
print("h = 2 * r * tan(theta)")