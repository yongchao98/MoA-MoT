import cmath

# Define the complex numbers for the given integrals.
# In Python, the imaginary unit is represented by 'j'.
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# According to the principle of deformation of contours, the integral over γ
# is the sum of the integrals over γ1 and γ2 because γ encircles both
# singularities z1 and z2 in the same direction (counter-clockwise).
integral_gamma = integral_gamma1 + integral_gamma2

# Extract integer parts for clean printing.
r1, i1 = int(integral_gamma1.real), int(integral_gamma1.imag)
r2, i2 = int(integral_gamma2.real), int(integral_gamma2.imag)
r_result, i_result = int(integral_gamma.real), int(integral_gamma.imag)

# Print the final equation showing the calculation.
print("The calculation is based on the formula: ∫γ f = ∫γ₁ f + ∫γ₂ f")
print("Substituting the given values:")
print(f"({r1} + {i1}i) + ({r2} + {i2}i) = {r_result} + {i_result}i")
