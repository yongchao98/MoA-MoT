import sympy

# Define symbolic variables for parametrization
u, v = sympy.symbols('u v')

# The problem is to calculate the flux of a vector field F through two yellow sides of a pyramid.
# Flux = Integral(F . dS)
# Due to symmetry, the flux through the two opposite yellow faces is identical.
# We will calculate the flux for one face and multiply the result by 2.

# --- Calculation for one yellow side (the face in the positive y direction) ---

# Parametrization of the triangular face with vertices at (0,0,4), (-1,1,0), and (1,1,0).
# r(u,v) = <u(1-v), 1-v, 4v> for u in [-1,1] and v in [0,1].
r = sympy.Matrix([u * (1 - v), 1 - v, 4 * v])

# Calculate the partial derivatives to find the surface element dS.
dr_du = sympy.diff(r, u)
dr_dv = sympy.diff(r, v)

# The normal vector dS is given by the cross product of the partial derivatives.
# The cross product dr_du x dr_dv results in an inward-pointing vector.
# We take its negative to get the required outward-pointing normal vector.
dS_vec = -dr_du.cross(dr_dv)

# Define the vector field F = <3x^3*y^2*z, 3x^2*y^3, z>.
x, y, z = sympy.symbols('x y z')
F_field = sympy.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

# Evaluate the vector field F on the parametrized surface.
x_s, y_s, z_s = r
F_on_surface = F_field.subs([(x, x_s), (y, y_s), (z, z_s)])

# The integrand is the dot product of F on the surface and the normal vector dS.
integrand = F_on_surface.dot(dS_vec)

# It's helpful to simplify the integrand before integration.
integrand_simplified = sympy.simplify(integrand)

# To find the flux, integrate the simplified expression over the bounds of u and v.
# First, integrate with respect to u from -1 to 1.
flux_after_u_integration = sympy.integrate(integrand_simplified, (u, -1, 1))

# Then, integrate the result with respect to v from 0 to 1.
flux_one_side = sympy.integrate(flux_after_u_integration, (v, 0, 1))

# The total flux is the sum of the fluxes through the two identical yellow sides.
total_flux = 2 * flux_one_side

# Print the final equation showing the flux for each side and the total.
print("The total energy flow is the sum of the flows through each of the two yellow sides.")
print("Final Equation:")
print(f"{flux_one_side} + {flux_one_side} = {total_flux}")