# Define the given integral values as complex numbers.
# In Python, the imaginary unit is represented by 'j'.
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# Based on the deformation principle, the integral over γ is the sum of the
# integrals over its two loops.
# The loop around z2 has the same orientation as γ2.
# The loop around z1 has the opposite orientation of γ1.
# So, ∫γ f = ∫γ₂ f - ∫γ₁ f.
result = integral_gamma2 - integral_gamma1

# Format the complex numbers for a clear output string matching mathematical notation.
def format_complex(c):
    # Use 'i' for the imaginary unit in the output string.
    return str(c).replace('j', 'i')

# Print the final equation with all the numbers.
print(f"The integral over γ is calculated as follows:")
print(f"∫γ f = ∫γ₂ f - ∫γ₁ f")
print(f"     = {format_complex(integral_gamma2)} - {format_complex(integral_gamma1)}")
print(f"     = {format_complex(result)}")
