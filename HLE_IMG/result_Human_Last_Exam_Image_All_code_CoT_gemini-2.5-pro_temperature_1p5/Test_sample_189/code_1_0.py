import cmath

# Given values for the integrals
integral_gamma1 = 3 + 4j
integral_gamma2 = 5 + 6j

# According to the deformation principle for contours,
# the integral over γ is the sum of the integrals over γ1 and γ2.
integral_gamma = integral_gamma1 + integral_gamma2

# Format the output to clearly show the calculation
# using the given numbers.
real1 = int(integral_gamma1.real)
imag1 = int(integral_gamma1.imag)
real2 = int(integral_gamma2.real)
imag2 = int(integral_gamma2.imag)

real_result = int(integral_gamma.real)
imag_result = int(integral_gamma.imag)

# Print the final equation
print(f"∫γ f = ∫γ1 f + ∫γ2 f")
print(f"∫γ f = ({real1} + {imag1}i) + ({real2} + {imag2}i)")
print(f"∫γ f = ({real1} + {real2}) + ({imag1} + {imag2})i")
print(f"∫γ f = {real_result} + {imag_result}i")
