import mpmath

# Set precision for the calculation
mpmath.mp.dps = 50

# Define the parameters for the Hurwitz zeta functions
s = 1.5
a_val = 1.5
b_val = 2/3
a_complex = mpmath.mpc(a_val, -b_val)

# Calculate the values of the Hurwitz zeta functions
zeta1 = mpmath.hurwitz(s, a_complex)
zeta2 = mpmath.hurwitz(s, mpmath.conj(a_complex))

# Calculate the final value of the integral
i = mpmath.mpc(0, 1)
sqrt_pi = mpmath.sqrt(mpmath.pi)
integral_value = i * sqrt_pi * (zeta1 + zeta2)

# Create the string for the formula
s_str = "3/2"
a_str = "3/2"
b_str = "2/3"
formula = f"i * sqrt(pi) * (zeta({s_str}, {a_str} - i*{b_str}) + zeta({s_str}, {a_str} + i*{b_str}))"

# Print the formula and the calculated value
print(f"The analytical value of the integral is given by the formula:")
print(formula)
print("\nWhere:")
print(f"i is the imaginary unit")
print(f"pi is the mathematical constant pi ({mpmath.pi})")
print(f"sqrt is the square root function")
print(f"zeta(s, a) is the Hurwitz zeta function")
print(f"zeta({s_str}, {a_str} - i*{b_str}) = {zeta1}")
print(f"zeta({s_str}, {a_str} + i*{b_str}) = {zeta2}")
print(f"\nThe numerical value is:")
print(f"{formula} = {integral_value}")

# Return the final numerical answer in the required format
final_answer = integral_value
# The final answer is a complex number, we format it as requested.
# As the real part is negligible (due to precision), we only output the imaginary part.
final_answer_str = f"{mpmath.nstr(final_answer.imag, 15)}j"
# <<<...>>> is not suitable for a complex number string output. 
# We print the value and put it in the required format if it was a single number.
# Let's provide the imaginary part as the final answer.
final_imag_part = final_answer.imag