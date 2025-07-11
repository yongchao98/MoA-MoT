import sympy

# Set up the variable for the polynomial
t = sympy.Symbol('t')

# The primitive degrees for the quaternionic reflection group associated with W(E_8)
degrees1 = [3, 5, 7]
degrees2 = [11, 13, 15]
degrees3 = [19, 21, 23]
degrees4 = [29, 31, 33]

# Construct the factors of the Poincare polynomial
P1 = 1
for d in degrees1:
    P1 *= (1 + t**d)

P2 = 1
for d in degrees2:
    P2 *= (1 + t**d)

P3 = 1
for d in degrees3:
    P3 *= (1 + t**d)
    
P4 = 1
for d in degrees4:
    P4 *= (1 + t**d)

# The full Poincare polynomial is the product of these factors
P = P1 * P2 * P3 * P4

# Expand the first part of the polynomial to show the calculation
P1_expanded = sympy.expand(P1)
print(f"The first factor of the polynomial is (1+t^{degrees1[0]})(1+t^{degrees1[1]})(1+t^{degrees1[2]}) = {P1_expanded}")

# The degree we are interested in is 9
k = 9

# Extract the coefficient of t^k
# As reasoned in the plan, any contribution to t^9 must come from the first factor,
# since all powers in other factors are > 9.
b_k = P1_expanded.coeff(t, k)

# The full equation is dim H^9(M, Q) = [t^9] P(t)
# We show the numbers involved in the definition of P(t) and the result.
print(f"The Poincare polynomial P(t) is the product of factors (1+t^d) for d in {degrees1 + degrees2 + degrees3 + degrees4}.")
print(f"The dimension of the ninth cohomology group is the coefficient of t^{k} in P(t).")
print(f"dim H^{k}(M, Q) = {b_k}")