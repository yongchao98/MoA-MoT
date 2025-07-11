import cmath

# Given values of the integrals
integral1 = 3 + 4j
integral2 = 5 + 6j

# The contour γ is homologous to the sum of contours γ₁ and γ₂.
# This means the integral over γ is the sum of the integrals over γ₁ and γ₂.
# ∫_γ f = ∫_γ₁ f + ∫_γ₂ f

# Calculate the sum of the two integrals
result_integral = integral1 + integral2

# Extract the real and imaginary parts for clear printing
r1, i1 = integral1.real, integral1.imag
r2, i2 = integral2.real, integral2.imag
r_res, i_res = result_integral.real, result_integral.imag

# Print the equation and the final answer
print("The integral over the contour γ is the sum of the integrals over contours γ₁ and γ₂.")
print("∫_γ f = ∫_γ₁ f + ∫_γ₂ f")
# Using int() to remove the .0 for cleaner output
print(f"∫_γ f = ({int(r1)} + {int(i1)}i) + ({int(r2)} + {int(i2)}i)")
print(f"∫_γ f = ({int(r1)} + {int(r2)}) + ({int(i1)} + {int(i2)})i")
print(f"∫_γ f = {int(r_res)} + {int(i_res)}i")
