import cmath

# Define the given values of the integrals as complex numbers.
# In Python, the imaginary unit is represented by 'j'.
integral_gamma1 = complex(3, 4)
integral_gamma2 = complex(5, 6)

# The contour γ winds counter-clockwise around z2 (positive direction)
# and clockwise around z1 (negative direction).
# By the principle of deformation of contours (a result of Cauchy's Integral Theorem),
# the integral over γ is the sum of the integrals over contours enclosing each singularity,
# weighted by the winding number.
#
# So, ∫_γ f(z) dz = (Winding number around z2) * ∫_γ2 f(z) dz + (Winding number around z1) * ∫_γ1 f(z) dz
# ∫_γ f(z) dz = (+1) * ∫_γ2 f(z) dz + (-1) * ∫_γ1 f(z) dz
# ∫_γ f(z) dz = ∫_γ2 f(z) dz - ∫_γ1 f(z) dz

# Calculate the integral over γ
integral_gamma = integral_gamma2 - integral_gamma1

# Extract real and imaginary parts for clear printing
ig1_real = integral_gamma1.real
ig1_imag = integral_gamma1.imag
ig2_real = integral_gamma2.real
ig2_imag = integral_gamma2.imag
ig_real = integral_gamma.real
ig_imag = integral_gamma.imag

# Print the calculation step-by-step
print("The integral over the contour γ is calculated by applying the deformation principle:")
print("∫γ f(z)dz = ∫γ2 f(z)dz - ∫γ1 f(z)dz")
print("\nSubstituting the given values:")
print(f"∫γ f(z)dz = ({int(ig2_real)} + {int(ig2_imag)}i) - ({int(ig1_real)} + {int(ig1_imag)}i)")
print("\nPerforming the subtraction of complex numbers:")
print(f"∫γ f(z)dz = ({int(ig2_real)} - {int(ig1_real)}) + ({int(ig2_imag)} - {int(ig1_imag)})i")
print("\nFinal Result:")
print(f"∫γ f(z)dz = {int(ig_real)} + {int(ig_imag)}i")