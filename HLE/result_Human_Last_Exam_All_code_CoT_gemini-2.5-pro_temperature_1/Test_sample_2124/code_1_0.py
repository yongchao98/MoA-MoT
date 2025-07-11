# The one-loop counter-term coefficients are proportional to a common factor,
# let's call it C = g^2 / (32 * pi^2 * epsilon). We can calculate the ratio
# using their relative values.

# From the fermion self-energy calculation:
# delta_Zx = -1 * C
dZx = -1.0

# From the fermion self-energy and mass renormalization:
# delta_Zmx = 3 * C
dZmx = 3.0

# From the Ward identity simplified by the condition delta_Z_phi = 0:
# delta_Zg = -1 * C
dZg = -1.0

# Calculate the ratio R = delta_Zx / (delta_Zg + delta_Zmx)
R = dZx / (dZg + dZmx)

# Output the numbers used in the final calculation
print(f"The relative value of the fermion field counter-term coefficient is: delta_Zx = {dZx}")
print(f"The relative value of the Yukawa coupling counter-term coefficient is: delta_Zg = {dZg}")
print(f"The relative value of the fermion mass counter-term coefficient is: delta_Zmx = {dZmx}")
print(f"The ratio R is calculated as: R = {dZx} / ({dZg} + {dZmx})")
print(f"The final result is: R = {R}")
<<< -0.5 >>>