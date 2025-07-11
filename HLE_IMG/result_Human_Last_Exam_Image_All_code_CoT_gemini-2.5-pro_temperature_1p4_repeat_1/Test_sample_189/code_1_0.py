import cmath

# Given integral values
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# The contour gamma is homotopic to gamma2 traversed counter-clockwise and
# gamma1 traversed clockwise.
# Thus, the integral over gamma is the integral over gamma2 minus the integral over gamma1.
result = integral_gamma2 - integral_gamma1

# Extract real and imaginary parts for clear printing
g1_real = int(integral_gamma1.real)
g1_imag = int(integral_gamma1.imag)
g2_real = int(integral_gamma2.real)
g2_imag = int(integral_gamma2.imag)
res_real = int(result.real)
res_imag = int(result.imag)

# Print the final calculation as an equation
print("The integral over the contour gamma is calculated as follows:")
print(f"∫γ f = ∫γ2 f - ∫γ1 f")
print(f"     = ({g2_real} + {g2_imag}i) - ({g1_real} + {g1_imag}i)")
print(f"     = ({g2_real} - {g1_real}) + ({g2_imag} - {g1_imag})i")
print(f"     = {res_real} + {res_imag}i")