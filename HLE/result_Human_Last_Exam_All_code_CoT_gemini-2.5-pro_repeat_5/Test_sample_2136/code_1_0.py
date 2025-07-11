# Coefficients from the PDE
# Viscosity term coefficient
nu = 3
# Value at the stationary point
# This is derived from the condition on the spatial gradient at (x0, tau)
# Given u_x(x0, tau) = -1, and for this type of solution u_x = -C^2.
# So, -1 = -C^2, which gives C = 1.
C = 1.0

# Formula for the integral for this specific type of solution
# I = (54/5) * nu * C^3
numerator = 54
denominator = 5
integral_value = (numerator / denominator) * nu * (C**3)

# Output the equation with the numbers plugged in
print(f"The integral is calculated using the formula: I = ({numerator}/{denominator}) * nu * C^3")
print(f"Plugging in the values nu = {nu} and C = {C}:")
print(f"I = ({numerator}/{denominator}) * {nu} * ({C})^3")
print(f"I = {integral_value}")
