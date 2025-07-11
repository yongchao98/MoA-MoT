# In this problem, we are asked to calculate the ratio of one-loop counter-terms
# in a Yukawa theory. The steps involve calculating the divergent parts of the
# fermion self-energy and the vertex correction.

# Let C be the common factor g^2 / (32 * pi^2 * epsilon).
# We express all counter-terms in units of C.

# From the fermion self-energy calculation, the counter-term for the fermion field is:
# delta_Z_x = g^2 / (32 * pi^2 * epsilon)
delta_Z_x_coeff = 1.0

# From the fermion self-energy calculation, the counter-term for the fermion mass is:
# delta_Z_m_x = -g^2 / (16 * pi^2 * epsilon) = -2 * C
delta_Z_m_x_coeff = -2.0

# From the vertex correction calculation, the counter-term for the Yukawa coupling is:
# delta_Z_g = -g^2 / (16 * pi^2 * epsilon) = -2 * C
delta_Z_g_coeff = -2.0

# The ratio R is defined as delta_Z_x / (delta_Z_g + delta_Z_m_x).
# The common factor C will cancel out.
numerator = delta_Z_x_coeff
denominator = delta_Z_g_coeff + delta_Z_m_x_coeff

# Calculate the final ratio R.
R = numerator / denominator

# We print the equation with the coefficients to show the calculation.
print(f"The ratio R is calculated as:")
print(f"R = (delta_Z_x) / (delta_Z_g + delta_Z_m_x)")
print(f"  = ({delta_Z_x_coeff} * C) / (({delta_Z_g_coeff} * C) + ({delta_Z_m_x_coeff} * C))")
print(f"  = {delta_Z_x_coeff} / ({delta_Z_g_coeff} + {delta_Z_m_x_coeff})")
print(f"  = {numerator} / {denominator}")
print(f"  = {R}")

print("\nFinal Answer:")
print(f"<<<{R}>>>")