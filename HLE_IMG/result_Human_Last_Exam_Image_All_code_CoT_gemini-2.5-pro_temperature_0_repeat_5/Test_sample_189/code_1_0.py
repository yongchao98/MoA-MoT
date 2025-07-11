# Define the values of the integrals as complex numbers
integral_g1 = complex(3, 4)
integral_g2 = complex(5, 6)

# The integral over gamma is the sum of the integrals over gamma1 and gamma2
integral_g = integral_g1 + integral_g2

# Helper function to format complex numbers for printing
def format_complex(c):
    real_part = int(c.real)
    imag_part = int(c.imag)
    if imag_part >= 0:
        return f"({real_part} + {imag_part}i)"
    else:
        return f"({real_part} - {abs(imag_part)}i)"

# Print the final equation with the numbers
print(f"By the deformation principle, the integral over γ is the sum of the integrals over γ1 and γ2.")
print(f"∫γ f = ∫γ1 f + ∫γ2 f")
print(f"∫γ f = {format_complex(integral_g1)} + {format_complex(integral_g2)}")
print(f"∫γ f = {int(integral_g.real)} + {int(integral_g.imag)}i")
