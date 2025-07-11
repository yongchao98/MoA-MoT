import cmath

# Given integral values
integral_gamma1 = complex(3, 4)
integral_gamma2 = complex(5, 6)

# The contour γ is equivalent to γ₂ (counter-clockwise around z₂) 
# and -γ₁ (clockwise around z₁).
# Therefore, the integral over γ is ∫f dγ = ∫f dγ₂ - ∫f dγ₁.
result = integral_gamma2 - integral_gamma1

# Custom function to format complex numbers to the "a + bi" notation
def format_complex(c):
    real = c.real
    imag = c.imag
    # Use integers since the inputs are integers
    if real.is_integer():
        real = int(real)
    if imag.is_integer():
        imag = int(imag)
    
    if imag == 0:
        return str(real)
    elif real == 0:
        return f"{imag}i"
    elif imag > 0:
        return f"{real} + {imag}i"
    else:
        # The minus sign is already part of the imaginary number
        return f"{real} - {-imag}i"

# Print the final equation with the numbers
print("The final calculation is:")
print(f"({format_complex(integral_gamma2)}) - ({format_complex(integral_gamma1)}) = {format_complex(result)}")
