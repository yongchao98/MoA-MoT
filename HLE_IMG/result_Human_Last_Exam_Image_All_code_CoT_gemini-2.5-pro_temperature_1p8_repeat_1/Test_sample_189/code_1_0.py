# Define the complex numbers representing the given integrals
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# The contour γ makes a counter-clockwise loop around z₁ and a counter-clockwise loop around z₂.
# By the principle of deformation of contours (or Cauchy's Integral Theorem for multiply-connected domains),
# the integral over γ is the sum of the integrals over contours γ₁ and γ₂, as γ is homologous to γ₁ + γ₂.
# ∫_γ f = ∫_γ₁ f + ∫_γ₂ f

# Perform the addition of the complex numbers
integral_gamma = integral_gamma1 + integral_gamma2

# Extract the real and imaginary parts of all numbers to construct the final output string
r1 = int(integral_gamma1.real)
i1 = int(integral_gamma1.imag)
r2 = int(integral_gamma2.real)
i2 = int(integral_gamma2.imag)
final_r = int(integral_gamma.real)
final_i = int(integral_gamma.imag)

# Print the final calculation as a full equation, including each number
print("Based on the principle of deformation of contours, the integral over γ is the sum of the integrals over γ₁ and γ₂.")
print(f"∫_γ f = ({r1} + {i1}i) + ({r2} + {i2}i)")
print(f"∫_γ f = ({r1} + {r2}) + ({i1} + {i2})i")
print(f"∫_γ f = {final_r} + {final_i}i")