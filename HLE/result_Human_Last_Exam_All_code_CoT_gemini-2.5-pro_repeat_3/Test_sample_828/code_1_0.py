import sympy

# 1. Define symbolic variables for the physical quantities and coordinates.
# 'a' and 'omega' are constants of the system.
# 'r' and 'theta' are the polar coordinates.
# 'T_cal' represents the script T symbol (mathcal{T}) from the problem.
a, omega, r, theta = sympy.symbols('a omega r theta', real=True, positive=True)
T_cal = sympy.Symbol('T_cal')

# 2. Define the given Cartesian tensor components T_ij.
# These components are expressed in terms of the polar angle theta.
T_xx = T_cal * (a * omega)**2 * sympy.sin(theta)**2
T_yy = T_cal * (a * omega)**2 * sympy.cos(theta)**2
T_xy = -T_cal * (a * omega)**2 * sympy.sin(theta) * sympy.cos(theta)

# 3. Calculate the partial derivatives for the coordinate transformation.
# We need the derivatives of Cartesian coordinates (x, y) with respect to theta.
# x = r*cos(theta), y = r*sin(theta)
dx_dtheta = sympy.diff(r * sympy.cos(theta), theta)
dy_dtheta = sympy.diff(r * sympy.sin(theta), theta)

# 4. Apply the tensor transformation law to find the T_thetatheta component in polar coordinates.
# T'_thetatheta = (dx/dtheta)^2 * T_xx + (dy/dtheta)^2 * T_yy + 2*(dx/dtheta)*(dy/dtheta)*T_xy
T_prime_thetatheta = (dx_dtheta**2 * T_xx +
                      dy_dtheta**2 * T_yy +
                      2 * dx_dtheta * dy_dtheta * T_xy)

# We simplify the resulting expression.
T_prime_thetatheta_simplified = sympy.simplify(T_prime_thetatheta)

# 5. The problem concerns a ring of radius 'a', so we evaluate the component on this ring by substituting r=a.
T_thetatheta_on_ring = T_prime_thetatheta_simplified.subs(r, a)

# 6. The problem asks to find K from the relation: T_thetatheta = a^2*sin(theta)^2*T_cal + K.
# We define the known part of this expression.
given_term_part1 = a**2 * sympy.sin(theta)**2 * T_cal

# 7. By equating our calculated result with the given form, we solve for K.
# K = T_thetatheta_on_ring - given_term_part1
K_solution = sympy.simplify(T_thetatheta_on_ring - given_term_part1)

# Print the final equation with all parts shown explicitly.
print("The transformed component T_thetatheta on the ring is:")
sympy.pprint(T_thetatheta_on_ring)
print("\nThis is equated to the expression from the problem: a**2*sin(theta)**2*T_cal + K")
print("\nSo, the full equation is:")
# Create a symbolic K for printing the equation
K_symbol = sympy.Symbol('K')
equation = sympy.Eq(T_thetatheta_on_ring, given_term_part1 + K_symbol)
sympy.pprint(equation)

print("\nSolving for K, we find:")
sympy.pprint(K_solution)