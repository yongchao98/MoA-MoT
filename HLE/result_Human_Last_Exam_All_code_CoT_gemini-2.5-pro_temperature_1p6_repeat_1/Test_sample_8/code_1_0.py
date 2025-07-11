import math

# Let alpha be a representative value in the interval (0, 1)
# The exact value does not change the outcome of the comparison.
alpha = 0.5

# The conormal exponent for the source function f is 2 + alpha.
# This corresponds to a decay of r^-(2+alpha) at infinity.
source_exponent = 2 + alpha

# The conormal exponent for the outgoing scattered wave is 1.
# This corresponds to a decay of r^-1 at infinity.
scattering_exponent = 1

# The resulting function's decay is determined by the slower of the two rates,
# which corresponds to the smaller of the two exponents.
final_exponent = min(source_exponent, scattering_exponent)

print(f"The problem is to find the conormal space for R(sigma)f, where f is in A^(2+alpha)(X).")
print(f"Let's use a sample value for alpha = {alpha}.")
print(f"The source term's conormal exponent is k_source = 2 + {alpha} = {source_exponent}.")
print(f"The scattered wave's conormal exponent is k_scatter = {scattering_exponent}.")
print("\nThe final conormal exponent, k_final, is the minimum of these two.")
print(f"The final equation is: k_final = min(k_source, k_scatter)")
print(f"Substituting the values: k_final = min({source_exponent}, {scattering_exponent})")
print(f"The result of this equation is: k_final = {final_exponent}")
print(f"\nTherefore, R(sigma)f belongs to the conormal space A^{final_exponent}(X).")
