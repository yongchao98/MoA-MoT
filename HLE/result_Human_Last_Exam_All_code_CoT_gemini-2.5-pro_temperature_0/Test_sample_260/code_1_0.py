from fractions import Fraction

# The orders of the cyclic groups G1, G2, G3
n1 = 5
n2 = 8
n3 = 2

# Calculate the Euler characteristic of G = G1 * G2 * G3
# chi(G) = chi(G1) + chi(G2) + chi(G3) - 2
chi_g1 = Fraction(1, n1)
chi_g2 = Fraction(1, n2)
chi_g3 = Fraction(1, n3)
chi_G = chi_g1 + chi_g2 + chi_g3 - 2

# Calculate the order of the homology group H = G1 + G2 + G3
order_H = n1 * n2 * n3

# The Euler characteristic of H is 1/|H|
chi_H = Fraction(1, order_H)

# Using the formula chi(G) = chi(K) * chi(H), where chi(K) = 1 - r
# We get 1 - r = chi(G) / chi(H) = chi(G) * |H|
one_minus_r = chi_G * order_H

# Solve for the rank r
r = 1 - one_minus_r

print("The rank 'r' of the kernel K is calculated using the formula r = 1 - (chi(G) * |H|).")
print("The numbers in the final equation are:")
print(f"chi(G) = 1/{n1} + 1/{n2} + 1/{n3} - 2 = {chi_G}")
print(f"|H| = {n1} * {n2} * {n3} = {order_H}")
print("\nSubstituting these values into the equation for r:")
print(f"r = 1 - ({chi_G} * {order_H})")
print(f"r = 1 - ({one_minus_r})")
# The result must be an integer, as rank is an integer.
print(f"The rank of the kernel is: {int(r)}")
