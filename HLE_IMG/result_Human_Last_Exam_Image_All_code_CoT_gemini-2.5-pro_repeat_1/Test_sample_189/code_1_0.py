import cmath

# Define the given values of the integrals as complex numbers.
# Python uses 'j' for the imaginary unit.
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# According to the deformation principle and the winding numbers of the contour gamma,
# the integral over gamma is calculated as:
# integral(gamma) = integral(gamma_2) - integral(gamma_1)
# This is because gamma winds +1 times around z2 and -1 times around z1.

# Perform the complex number subtraction.
result = integral_gamma2 - integral_gamma1

# Print the equation showing all the numbers involved in the calculation.
print("The calculation for the integral over gamma is:")
# Note: Python's complex number representation is (real+imagj)
print(f"({integral_gamma2}) - ({integral_gamma1}) = {result}")

# To make the output clearer with 'i', we can format it.
def format_complex(c):
    # Format to use 'i' and handle signs nicely
    real_part = c.real
    imag_part = c.imag
    if imag_part >= 0:
        return f"{real_part} + {imag_part}i"
    else:
        return f"{real_part} - {-imag_part}i"

print("\nIn standard mathematical notation (using 'i'):")
print(f"({format_complex(integral_gamma2)}) - ({format_complex(integral_gamma1)}) = {format_complex(result)}")
