# The separatrix is described by the algebraic equation d = -u**2.
# We can rewrite this in the standard form: d + u**2 = 0.

# Define the numerical components of the equation in the form:
# (coeff_d) * d + (coeff_u) * u**(power_u) = constant
coeff_d = 1
coeff_u = 1
power_u = 2
constant = 0

# Print the final equation with each number clearly shown.
print("The equation for the separatrix is:")
print(f"({coeff_d})*d + ({coeff_u})*u**({power_u}) = {constant}")
