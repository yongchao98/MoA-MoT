import math

# Step 1: Define the connective constant of the graph G.
# The graph G is constructed from an infinite ladder graph with additional diagonal edges.
# This structure is known in statistical mechanics as the martini lattice.
# The connective constant (mu) of a graph is the exponential growth rate
# of the number of self-avoiding walks. For the martini lattice, this value is
# known to be mu = sqrt(2 + sqrt(2)).

mu = math.sqrt(2 + math.sqrt(2))
print(f"The connective constant mu is sqrt(2 + sqrt(2)) ≈ {mu:.8f}")

# Step 2: Derive the minimal polynomial for mu.
# Let x = mu. The derivation proceeds as follows:
# x = sqrt(2 + sqrt(2))
# Square both sides:
# x^2 = 2 + sqrt(2)
# Isolate the remaining square root:
# x^2 - 2 = sqrt(2)
# Square both sides again to eliminate the root:
# (x^2 - 2)^2 = 2
# Expand the equation:
# x^4 - 4*x^2 + 4 = 2
# Rearrange to the standard polynomial form P(x) = 0:
# x^4 - 4*x^2 + 2 = 0
# This gives us a candidate polynomial P(x) = x^4 - 4*x^2 + 2.

# Step 3: Verify that mu is a root of this polynomial.
# We can plug the numerical value of mu into the polynomial expression.
# The result should be very close to zero, with any small difference due to floating-point errors.
poly_val = mu**4 - 4 * mu**2 + 2
print(f"Verifying the polynomial: mu^4 - 4*mu^2 + 2 = {poly_val:.8f}")
print("The value is effectively zero, confirming mu is a root.")

# Step 4: Justify why this is the minimal polynomial.
# For P(x) = x^4 - 4*x^2 + 2 to be the minimal polynomial, it must be irreducible
# over the field of rational numbers (Q). We can use Eisenstein's criterion
# with the prime number p = 2.
# The coefficients are: a_4=1, a_3=0, a_2=-4, a_1=0, a_0=2.
# 1. p=2 divides all coefficients except the leading one (a_3, a_2, a_1, a_0). This is true.
# 2. p=2 does not divide the leading coefficient a_4=1. This is true.
# 3. p^2=4 does not divide the constant term a_0=2. This is true.
# Since all conditions of Eisenstein's criterion are met, the polynomial is irreducible over Q.
# Therefore, x^4 - 4*x^2 + 2 is the minimal polynomial.

# Step 5: Output the final answer in the required format.
print("\nThe minimal polynomial for the connective constant is: x^4 - 4*x^2 + 2 = 0")
print("\nThe numbers (non-zero coefficients) in the final equation are:")
coeff_x4 = 1
coeff_x2 = -4
coeff_x0 = 2
print(coeff_x4)
print(coeff_x2)
print(coeff_x0)