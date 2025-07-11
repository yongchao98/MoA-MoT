import mpmath

# Set the desired precision
mpmath.mp.dps = 50

# Define the constants from the derived formula for J
# J = 2 * sqrt(pi) * Re(zeta(s, q))
two = mpmath.mpf(2)
three = mpmath.mpf(3)
s = three / two
q_real = three / two
q_imag = -two / three
q = mpmath.mpc(q_real, q_imag)

# Calculate the components of the formula
sqrt_pi = mpmath.sqrt(mpmath.pi)
hurwitz_zeta_val = mpmath.zeta(s, q)
real_part_zeta = mpmath.re(hurwitz_zeta_val)

# Calculate the value of J
J = two * sqrt_pi * real_part_zeta

# The final answer for the integral I is i*J
# To satisfy the output format requirement, we will print the numbers in the final expression for J.
# J = 2 * sqrt(pi) * Re(zeta(1.5, 1.5 - 2/3i))
print("The analytical formula for the real part J is:")
print(f"J = {two} * sqrt(pi) * Re(zeta({s}, {q_real} + ({q_imag})*i))")
print("\nCalculating the numerical value of J:")
print(f"J = {two} * {sqrt_pi} * {real_part_zeta}")
print(f"J_value = {J}")

print("\nThe value of the original integral is I = i * J.")
print(f"Therefore, the final analytical value is approximately {J}*i.")
