# Define the complex numbers representing the given integrals.
integral_g1 = 3 + 4j
integral_g2 = 5 + 6j

# According to the deformation principle in complex analysis,
# the integral over the contour γ is the sum of the integrals over γ1 and γ2,
# as γ is homologous to γ1 + γ2.
integral_g = integral_g1 + integral_g2

# Extract the real and imaginary parts to display the calculation clearly.
# We cast to int because the inputs are integers.
g1_real = int(integral_g1.real)
g1_imag = int(integral_g1.imag)
g2_real = int(integral_g2.real)
g2_imag = int(integral_g2.imag)
g_real = int(integral_g.real)
g_imag = int(integral_g.imag)

# Print the final equation with all the numbers.
print(f"Let I(c) be the integral of f over a contour c.")
print(f"Based on the principle of deformation of contours, we have:")
print(f"I(γ) = I(γ1) + I(γ2)")
print(f"I(γ) = ({g1_real} + {g1_imag}i) + ({g2_real} + {g2_imag}i)")
print(f"I(γ) = ({g1_real} + {g2_real}) + ({g1_imag} + {g2_imag})i")
print(f"I(γ) = {g_real} + {g_imag}i")