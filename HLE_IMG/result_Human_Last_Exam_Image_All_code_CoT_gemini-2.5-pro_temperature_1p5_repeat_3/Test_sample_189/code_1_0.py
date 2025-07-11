import cmath

# Given values for the integrals
integral_g1 = complex(3, 4)
integral_g2 = complex(5, 6)

# The integral over the contour gamma is the integral over gamma_1 minus the integral over gamma_2
# due to the orientation of the loops in gamma.
# The loop in gamma around z1 is counter-clockwise, same as gamma_1.
# The loop in gamma around z2 is clockwise, opposite to gamma_2.
integral_g = integral_g1 - integral_g2

# We need to print the final equation including all numbers.
# To format complex numbers nicely (e.g., 3 + 4i instead of (3+4j))
def format_complex(c):
    real = c.real
    imag = c.imag
    # Handle integer cases to avoid .0
    real_part = f"{int(real)}" if real == int(real) else f"{real}"
    imag_part = f"{int(imag)}" if imag == int(imag) else f"{imag}"

    if imag == 0:
        return f"{real_part}"
    elif real == 0:
        return f"{imag_part}i"
    elif imag < 0:
        return f"{real_part} - {abs(imag_part)}i"
    else:
        return f"{real_part} + {imag_part}i"

# Print the final calculation and result
print(f"The integral over γ can be calculated by deforming the contour.")
print(f"∫_γ f = ∫_γ₁ f - ∫_γ₂ f")
print(f"∫_γ f = ({format_complex(integral_g1)}) - ({format_complex(integral_g2)})")
print(f"∫_γ f = {format_complex(integral_g)}")

# The final result in a simplified format for the final answer block.
# We just need the final value, not the whole equation.
final_answer = format_complex(integral_g)

# The line below is not needed for the user but useful for verification
# print(f"Final answer: {final_answer}")