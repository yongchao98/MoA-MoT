import cmath

# The values of the integrals are given as complex numbers.
# We use Python's complex type, where 'j' represents the imaginary unit.
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# From the analysis of the contours, the integral over γ is the integral over γ2
# minus the integral over γ1. This is because γ loops counter-clockwise around z2
# (same as γ2) and clockwise around z1 (opposite of γ1).
# ∫γ f = ∫γ2 f - ∫γ1 f
result = integral_gamma2 - integral_gamma1

# The following lines format the complex numbers for a clear output.
# The format specifier "{:.0f}" removes decimal points for integers.
# The format specifier "{:+.0f}" ensures a sign is always shown for the imaginary part.
str_gamma1 = f"({integral_gamma1.real:.0f}{integral_gamma1.imag:+.0f}i)"
str_gamma2 = f"({integral_gamma2.real:.0f}{integral_gamma2.imag:+.0f}i)"
str_result = f"{result.real:.0f}{result.imag:+.0f}i"

# Print the final equation with all the numbers.
print("The calculation is based on the decomposition of the contour γ.")
print(f"∫γ f = ∫γ2 f - ∫γ1 f")
print(f"∫γ f = {str_gamma2} - {str_gamma1}")
print(f"∫γ f = {str_result}")