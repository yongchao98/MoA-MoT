import math

# Parameters from the omitted simulation
m = 1.5
E = 2.0
delta_z = 1.0
n0 = 9

# The condition for the calculation
V = 2 * E

# Check if the formula condition V > 2m holds
if not V > 2 * m:
    print("Warning: Formula condition V > 2m is not met.")

# Calculate the momentum inside the barrier
# q = sqrt((E - V)^2 - m^2) which simplifies to sqrt(E^2 - m^2) for V=2E
q = math.sqrt(E**2 - m**2)

# Calculate 1/|t|^2 = 1 + sin^2(q * delta_z)
inv_t_squared = 1 + (math.sin(q * delta_z))**2

# Calculate the final result
result = n0 * inv_t_squared

# Print the equation with all the numbers
print(f"Result = n0 / |t|^2")
print(f"n0 = {n0}")
print(f"1/|t|^2 = 1 + sin^2(sqrt(E^2 - m^2) * delta_z)")
print(f"1/|t|^2 = 1 + sin^2(sqrt({E}^2 - {m}^2) * {delta_z})")
print(f"1/|t|^2 = 1 + sin^2(sqrt({E**2} - {m**2}))")
print(f"1/|t|^2 = 1 + sin^2(sqrt({E**2 - m**2}))")
print(f"1/|t|^2 = {inv_t_squared}")
print(f"Result = {n0} * {inv_t_squared}")
print(f"Result = {result}")
