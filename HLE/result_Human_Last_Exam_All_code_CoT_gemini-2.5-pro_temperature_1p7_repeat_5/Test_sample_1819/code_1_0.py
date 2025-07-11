import sympy

# Define symbolic variables
x, y = sympy.symbols('x y')

# State the assumption about which faces are yellow.
print("This solution calculates the energy flow through the two yellow sides of the square pyramid.")
print("The pyramid has four side faces. 'Interspersed' coloring means opposite faces have the same color.")
print("This leads to two possibilities for the pair of yellow faces.")
print("We will assume the yellow faces are the 'front' (in the y < 0 region) and 'back' (in the y > 0 region) sides.")
print("-" * 30)

# --- Calculation for the back face ---
print("1. Flux through the back face (where y > 0):")
# The integrand is F · dS, which we found to be 12*x**2*y**3 + 4 - 4*y
integrand_back = 12*x**2*y**3 + 4 - 4*y
print(f"   - The integrand is F · dS = {integrand_back}")

# Integrate with respect to x
inner_integral_back = sympy.integrate(integrand_back, (x, -y, y))
inner_integral_back_simplified = sympy.simplify(inner_integral_back)
print(f"   - Integrating with respect to x from -y to y gives: {inner_integral_back_simplified}")

# Integrate with respect to y
flux_back = sympy.integrate(inner_integral_back_simplified, (y, 0, 1))
print(f"   - Integrating the result with respect to y from 0 to 1 gives the flux for the back face:")
print(f"   - Flux_back = {flux_back}")
print("-" * 30)


# --- Calculation for the front face ---
print("2. Flux through the front face (where y < 0):")
# The integrand is F · dS, which is -12*x**2*y**3 + 4 + 4*y
integrand_front = -12*x**2*y**3 + 4 + 4*y
print(f"   - The integrand is F · dS = {integrand_front}")

# Integrate with respect to x
inner_integral_front = sympy.integrate(integrand_front, (x, y, -y))
inner_integral_front_simplified = sympy.simplify(inner_integral_front)
print(f"   - Integrating with respect to x from y to -y gives: {inner_integral_front_simplified}")

# Integrate with respect to y
flux_front = sympy.integrate(inner_integral_front_simplified, (y, -1, 0))
print(f"   - Integrating the result with respect to y from -1 to 0 gives the flux for the front face:")
print(f"   - Flux_front = {flux_front}")
print("-" * 30)


# --- Total Flux Calculation ---
total_flux = flux_back + flux_front
print("3. Total flux through both yellow sides:")
print("The total energy flow is the sum of the fluxes through the two faces.")
print(f"Total Flow = {flux_back} + {flux_front} = {total_flux}")
