# Define the complex numbers for the given integrals.
# In Python, the imaginary unit is represented by 'j'.
integral_1 = 3 + 4j
integral_2 = 5 + 6j

# According to the principle of deformation of contours, the integral over γ
# is the sum of the integrals over γ1 and γ2.
integral_gamma = integral_1 + integral_2

# Extract the real and imaginary parts for formatted printing.
real1, imag1 = int(integral_1.real), int(integral_1.imag)
real2, imag2 = int(integral_2.real), int(integral_2.imag)
real_res, imag_res = int(integral_gamma.real), int(integral_gamma.imag)

# Print the final equation with all the numbers.
print(f"The calculation is based on the formula: ∫_γ f = ∫_γ1 f + ∫_γ2 f")
print(f"Plugging in the given values:")
print(f"∫_γ f = ({real1} + {imag1}i) + ({real2} + {imag2}i)")
print(f"∫_γ f = ({real1} + {real2}) + ({imag1} + {imag2})i")
print(f"∫_γ f = {real_res} + {imag_res}i")