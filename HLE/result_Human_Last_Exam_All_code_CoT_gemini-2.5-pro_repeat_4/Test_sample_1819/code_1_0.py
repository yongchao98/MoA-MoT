import sympy

# Define symbols
x, y, z = sympy.symbols('x y z')

# Define the vector field F
F_vec = sympy.Matrix([3 * x**3 * y**2 * z, 3 * x**2 * y**3, z])

print("This script calculates the energy flow (flux) through two opposite sides of a pyramid.")
print("We assume the 'yellow sides' are the faces with their base on the y=1 and y=-1 lines.")
print("-" * 50)

# --- Face 1 (S1, right side, y > 0) ---
# The surface is defined by z = 4 - 4y.
# The outward normal vector for the surface element dS is (0, 4, 1) dx dy.
z_s1 = 4 - 4*y
N1 = sympy.Matrix([0, 4, 1])

# The dot product F . N1, with z substituted
integrand1 = F_vec.dot(N1).subs(z, z_s1)

# The integral is over the projection on the xy-plane: a triangle
# with limits y from 0 to 1, and for each y, x from -y to y.
flux1 = sympy.integrate(integrand1, (x, -y, y), (y, 0, 1))

print("Calculation for the first yellow face (y > 0):")
print(f"The integrand is: {integrand1}")
print(f"The flux is: {flux1}")
print("-" * 50)

# --- Face 2 (S2, left side, y < 0) ---
# The surface is defined by z = 4 + 4y.
# The outward normal vector for the surface element dS is (0, -4, 1) dx dy.
z_s2 = 4 + 4*y
N2 = sympy.Matrix([0, -4, 1])

# The dot product F . N2, with z substituted
integrand2 = F_vec.dot(N2).subs(z, z_s2)

# The integral is over the projection on the xy-plane: a triangle
# with limits y from -1 to 0, and for each y, x from y to -y.
flux2 = sympy.integrate(integrand2, (x, y, -y), (y, -1, 0))

print("Calculation for the second yellow face (y < 0):")
print(f"The integrand is: {integrand2}")
print(f"The flux is: {flux2}")
print("-" * 50)

# --- Total Flow ---
total_flux = flux1 + flux2
print("The total energy flow is the sum of the fluxes through the two yellow faces.")
print(f"Final Equation: {flux1} + {flux2} = {total_flux}")