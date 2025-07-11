# Define the complex numbers for the given integrals
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# The integral over γ is the sum of the integrals over γ1 and γ2
# due to the principle of deformation of contours.
integral_gamma = integral_gamma1 + integral_gamma2

# Helper function to format complex numbers for printing
def format_complex(c):
    # Use .real and .imag attributes for the parts of the complex number
    # Format to remove trailing .0 for integers
    real_part = f"{c.real:.0f}"
    imag_part = f"{abs(c.imag):.0f}"
    
    if c.imag >= 0:
        return f"({real_part} + {imag_part}i)"
    else:
        return f"({real_part} - {imag_part}i)"

# Print the final equation showing all the numbers
print(f"∫γ f = ∫γ1 f + ∫γ2 f")
print(f"∫γ f = {format_complex(integral_gamma1)} + {format_complex(integral_gamma2)}")
print(f"∫γ f = {format_complex(integral_gamma)}")