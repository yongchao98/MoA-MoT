import scipy.constants

# 1. Define the given values and constants.
# Uncertainty in position (Δx) in meters. Given as 10 pm.
delta_x_pm = 10
delta_x_m = delta_x_pm * 1e-12

# Bohr radius (a₀), the radius of the first Bohr orbit, from scipy constants.
a_0 = scipy.constants.physical_constants["Bohr radius"][0]

# 2. The ratio Δp/p simplifies to a₀ / (2 * Δx).
# Δp = ħ / (2 * Δx) from Heisenberg's Uncertainty Principle.
# p = ħ / a₀ for the first Bohr orbit.
# Ratio = Δp / p = (ħ / (2 * Δx)) / (ħ / a₀) = a₀ / (2 * Δx).

# 3. Calculate the ratio.
ratio = a_0 / (2 * delta_x_m)

# 4. Print the final equation and the result.
print(f"The ratio is calculated as: Bohr radius / (2 * uncertainty in position)")
print(f"Ratio = {a_0} m / (2 * {delta_x_m} m)")
print(f"Ratio = {ratio}")
