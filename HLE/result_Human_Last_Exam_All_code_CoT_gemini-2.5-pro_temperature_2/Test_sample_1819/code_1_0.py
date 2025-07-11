import sympy

# Define the symbols for our variables
x, y = sympy.symbols('x y')

# --- Flux through the Front Face (S_front) ---

# The front face is on the plane z = 4 - 4y.
z_front = 4 - 4*y

# The vector field F is (3x^3y^2z, 3x^2y^3, z).
# The outward normal vector n is (0, 4, 1).
# We calculate the dot product F . n and substitute z.
# F . n = 12*x**2*y**3 + z
integrand_front = 12*x**2*y**3 + z_front

# The projection of the front face on the xy-plane is a triangle
# with vertices (0,0), (1,1), (-1,1).
# We integrate over this region: y from 0 to 1, and x from -y to y.
integral_in_x_front = sympy.integrate(integrand_front, (x, -y, y))
flux_front = sympy.integrate(integral_in_x_front, (y, 0, 1))


# --- Flux through the Back Face (S_back) ---

# The back face is on the plane z = 4 + 4y.
z_back = 4 + 4*y

# The outward normal vector n is (0, -4, 1).
# F . n = -12*x**2*y**3 + z
integrand_back = -12*x**2*y**3 + z_back

# The projection of the back face on the xy-plane is a triangle
# with vertices (0,0), (-1,-1), (1,-1).
# We integrate over this region: y from -1 to 0, and x from y to -y.
integral_in_x_back = sympy.integrate(integrand_back, (x, y, -y))
flux_back = sympy.integrate(integral_in_x_back, (y, -1, 0))

# --- Total Flux ---
# The total flux is the sum of the fluxes through the two yellow faces.
total_flux = flux_front + flux_back

print("Assuming the 'front' and 'back' faces are yellow:")
print(f"Flux through front face: {flux_front}")
print(f"Flux through back face: {flux_back}")
print("The total energy flow is the sum of these fluxes.")
# The final equation output
print(f"{flux_front} + {flux_back} = {total_flux}")
