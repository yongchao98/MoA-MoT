import sympy

# Define the symbols for integration
u, v = sympy.symbols('u v')

# The flux calculation for the two yellow sides simplifies to twice the value
# of a single integral due to symmetry. We calculate this integral in two parts.

# Part 1 of the integrand
integrand1 = 24 * (v - u)**2 * (u + v)**3

# Part 2 of the integrand
integrand2 = 8 - 8 * (u + v)

# Calculate the first integral (I1)
# The integration limits are v from 0 to 1-u, and u from 0 to 1.
integral1_val = sympy.integrate(integrand1, (v, 0, 1 - u), (u, 0, 1))

# Calculate the second integral (I2)
integral2_val = sympy.integrate(integrand2, (v, 0, 1 - u), (u, 0, 1))

# The flux through one side is the sum of the two parts
flux_one_side = integral1_val + integral2_val

# The total flux is twice the flux through one side
total_flux = 2 * flux_one_side

print("The calculation of the total energy flow is broken down as follows:")
print(f"The flux is given by 2 * (I1 + I2), where I1 and I2 are two separate integrals.")
print("-" * 30)
# Remember in the final code you still need to output each number in the final equation!
print(f"Value of Integral I1 = {integral1_val}")
print(f"Value of Integral I2 = {integral2_val}")
print("-" * 30)
print("The final equation for the total flow is:")
print(f"Total Flow = 2 * ({integral1_val} + {integral2_val})")
print(f"Total Flow = 2 * ({flux_one_side})")
print(f"Total Flow = {total_flux}")
print(f"The numerical value is approximately: {float(total_flux):.4f}")
