import sympy

# Define a common factor for the counter-terms.
# Let C = g^2 / (32 * pi^2 * epsilon)
# We can set this common factor to 1 for the calculation of the ratio.
C = 1

# From the one-loop calculation of the fermion self-energy in the MS-bar scheme, we get:
delta_Zx = 1 * C
print(f"The coefficient for the fermion field counter-term is: delta_Zx = {delta_Zx} * C")

# From the same calculation, we get the coefficient for the fermion mass counter-term:
delta_Zmx = -2 * C
print(f"The coefficient for the fermion mass counter-term is: delta_Zmx = {delta_Zmx} * C")

# From the vertex correction and the condition delta_Z_phi = 0, a Ward-like identity gives:
delta_Zg = 1 * C
print(f"The coefficient for the Yukawa coupling counter-term is: delta_Zg = {delta_Zg} * C")
print("-" * 20)

# Now, we calculate the ratio R
# R = delta_Zx / (delta_Zg + delta_Zmx)

# Calculate the denominator
denominator = delta_Zg + delta_Zmx

# Calculate the ratio R
R = delta_Zx / denominator

# Print the final equation and the result
print(f"The ratio R is calculated as:")
print(f"R = delta_Zx / (delta_Zg + delta_Zmx)")
print(f"R = {delta_Zx} / ({delta_Zg} + ({delta_Zmx}))")
print(f"R = {delta_Zx} / ({denominator})")
print(f"R = {R}")
