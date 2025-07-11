import mpmath

# Set the precision for the calculation
mpmath.mp.dps = 30

# Define the constants and parameters for the expression
# I = 2i * sqrt(pi) * Re(zeta(s, q))
s = mpmath.mpf('1.5')
q_real = mpmath.mpf('1.5')
q_imag = mpmath.mpf('-2/3')
q = mpmath.mpc(q_real, q_imag)

# Calculate the components of the formula
sqrt_pi = mpmath.sqrt(mpmath.pi)
zeta_val = mpmath.zeta(s, q)
real_part_zeta = zeta_val.real

# Calculate the final value of the integral
# The integral is purely imaginary
imaginary_part = 2 * sqrt_pi * real_part_zeta
final_value = mpmath.j * imaginary_part

# Print the equation with the computed values
print("The analytical formula for the integral is I = 2i * sqrt(pi) * Re(zeta(s, q))")
print(f"With s = {s}, q = {q}")
print("\nCalculating the components:")
print(f"sqrt(pi) approx {sqrt_pi}")
print(f"Re(zeta({s}, {q})) approx {real_part_zeta}")
print("\nThe final equation with numerical values is:")
print(f"I = 2j * {sqrt_pi} * {real_part_zeta}")
print(f"I approx {final_value}")

# The final answer is the imaginary part of the complex number
final_answer_val = final_value.imag