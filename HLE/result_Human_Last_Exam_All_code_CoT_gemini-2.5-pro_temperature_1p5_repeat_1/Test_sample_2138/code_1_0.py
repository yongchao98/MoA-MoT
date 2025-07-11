import mpmath

# Set the desired precision
mpmath.mp.dps = 50

# Define the parameters for the Hurwitz Zeta function
s = 1.5
# q = 3/2 + 2i/3
real_part_q = 1.5
imag_part_q = 2/3
q = mpmath.mpc(real_part_q, imag_part_q)

# Calculate the Hurwitz Zeta function zeta(s, q)
zeta_val = mpmath.hurwitz(s, q)

# The value of the integral is I = i * sqrt(pi) * (zeta(s, q) + zeta(s, q_conjugate))
# This simplifies to I = i * sqrt(pi) * 2 * Re(zeta(s, q))
# where Re is the real part.
sqrt_pi = mpmath.sqrt(mpmath.pi)
i = mpmath.j

integral_value = i * sqrt_pi * 2 * zeta_val.real

# Print the parameters used in the final expression
print(f"The integral is evaluated using the expression: 2 * i * sqrt(pi) * Re(zeta(s, q))")
print(f"s = {s}")
print(f"q = {real_part_q} + {imag_part_q}*i")

# Print the final result
print("\nThe analytical value of the integral is:")
print(integral_value)