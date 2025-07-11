# Define the given integral values as complex numbers
integral_gamma1 = complex(3, 4)
integral_gamma2 = complex(5, 6)

# The integral over γ is the sum of the integrals over γ1 and γ2
integral_gamma = integral_gamma1 + integral_gamma2

# Extract real and imaginary parts for clear printing
r1, i1 = integral_gamma1.real, integral_gamma1.imag
r2, i2 = integral_gamma2.real, integral_gamma2.imag
r_final, i_final = integral_gamma.real, integral_gamma.imag

# Print the final equation with all numbers
print(f"Based on the principle of deformation of contours, the integral over γ is the sum of the integrals over γ1 and γ2.")
print(f"∫(γ) f = ∫(γ1) f + ∫(γ2) f")
print(f"∫(γ) f = ({int(r1)} + {int(i1)}i) + ({int(r2)} + {int(i2)}i)")
print(f"∫(γ) f = ({int(r1)} + {int(r2)}) + ({int(i1)} + {int(i2)})i")
print(f"∫(γ) f = {int(r_final)} + {int(i_final)}i")