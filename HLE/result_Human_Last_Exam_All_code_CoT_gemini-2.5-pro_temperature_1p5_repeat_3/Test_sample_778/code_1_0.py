import math

# Define the parameters of the problem as provided.
n = 54669191328000
m = 1470432000
a = 991
b = 223

# Based on the theory of inserters of functors between delooped cyclic groups,
# the groupoid cardinality is given by the simple fraction m / n.
#
# A brief outline of the derivation:
# 1. The objects of the inserter groupoid are the 'm' elements of Z/mZ.
# 2. These objects are partitioned into d = gcd(b-a, m) connected components.
# 3. The size of the automorphism group for any object is |Aut| = (n * d) / m. This is guaranteed
#    to be an integer by the problem's premise that F and G are well-defined homomorphisms.
# 4. The groupoid cardinality is the number of components divided by the automorphism group size:
#    Cardinality = d / |Aut| = d / ((n * d) / m) = m / n.

# We now calculate this value and simplify the resulting fraction.

# To simplify the fraction m / n, we find the greatest common divisor (GCD).
greatest_common_divisor = math.gcd(m, n)

# The simplified numerator and denominator are obtained by dividing m and n by their GCD.
m_simplified = m // greatest_common_divisor
n_simplified = n // greatest_common_divisor

# Finally, we output the result in the requested format, showing the full equation.
print(f"The groupoid cardinality is the fraction m/n.")
print(f"The equation with the given numbers is:")
print(f"{m} / {n} = {m_simplified} / {n_simplified}")

# For completeness, we can also show the decimal value.
# final_value = m / n
# print(f"As a decimal, this is approximately: {final_value}")
