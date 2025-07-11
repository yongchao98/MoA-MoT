import sympy

# The problem of finding the maximum value of c_3 has been reduced to
# maximizing the function g(z1, z2) over the domain z1 in [-1, 0] and z2 in [0, 1].
# This reduction is based on the theory of moment spaces and extreme points of convex sets of measures.

# Define the symbolic variables
z1, z2 = sympy.symbols('z1 z2')

# Define the function to maximize, as derived in the thinking process
g = -z1 * z2 * (z1 + z2)

# The maximum of a function on a closed set occurs either at a critical point in the interior
# or on the boundary.
# A check for interior critical points (by setting partial derivatives to zero) reveals
# that the only critical point is at (0,0), which is on the boundary.
# So, we only need to check the boundary of the domain [-1, 0] x [0, 1].

# On boundaries z1=0 or z2=0, the function g is 0.

# Consider boundary where z2 = 1 and z1 is in [-1, 0]:
g_boundary = g.subs(z2, 1) # This gives -z1*(z1+1) = -z1**2 - z1

# To find the maximum of this function of z1 on [-1, 0], we take the derivative.
dg_boundary_dz1 = sympy.diff(g_boundary, z1)

# Solve for the critical point where the derivative is zero.
critical_points_z1 = sympy.solve(dg_boundary_dz1, z1)
opt_z1 = critical_points_z1[0]

# The maximum value of g is attained at this critical point (z1=-1/2, z2=1),
# as the function value is 0 at the endpoints of this boundary segment (z1=-1, z1=0).
max_g_value = g_boundary.subs(z1, opt_z1)

# The maximum value of c_3 is 35/2 times the maximum value of g.
c3_factor = sympy.Rational(35, 2)
max_c3 = c3_factor * max_g_value

# Print the results, including the final equation
print("The maximization problem for c_3 reduces to maximizing the function:")
print(f"g(z1, z2) = {g}")
print("over the domain z1 in [-1, 0] and z2 in [0, 1].\n")
print(f"The maximum value of g(z1, z2) is found to be {max_g_value}, occurring at z1={opt_z1}, z2=1.")
print("The maximum value of c_3 is then calculated as:")
print(f"max_c3 = (35/2) * max_g_value")
# Show the numbers in the final equation
print(f"max_c3 = {c3_factor} * {max_g_value} = {max_c3}")
