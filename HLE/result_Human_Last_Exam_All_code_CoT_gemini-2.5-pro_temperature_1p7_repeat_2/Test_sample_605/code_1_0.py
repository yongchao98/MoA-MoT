import math
from fractions import Fraction
from itertools import combinations

# 1. Define the weights from the problem statement.
weights = [22, 29, 49, 50, 75]

# 2. Calculate the degree 'd' for the Calabi-Yau manifold. It is the sum of the weights.
# We assume d = sum(weights) based on the Calabi-Yau condition, acknowledging the
# likely typo in the provided polynomial.
d = sum(weights)

# 3. Calculate the product of weights, P, and the sum of pairwise products, S2.
P = 1
for w in weights:
    P *= w

S2 = 0
for i, j in combinations(weights, 2):
    S2 += i * j

# 4. Calculate the second Chern number c2.H using the formula S2 * (d / P).
# We use Python's Fraction for exact arithmetic.
c2_H = Fraction(S2 * d, P)

# 5. Calculate the Crawley-Nordström invariant, which is (1/2) * c2.H.
CN_invariant = Fraction(1, 2) * c2_H

# 6. Print the components of the calculation and the final result as requested.
print("The Crawley-Nordström invariant `c` is calculated using the formula:")
print("c = (1/2) * (sum_{i<j} w_i*w_j) * (sum_k w_k) / (product_l w_l)")
print("")
print("For the given weights w = (22, 29, 49, 50, 75):")
print(f"The degree d = sum_k w_k = {d}")
print(f"The sum of products S2 = sum_{{i<j}} w_i*w_j = {S2}")
print(f"The product of weights P = product_l w_l = {P}")
print("")
print("The second Chern number c2.H = S2 * d / P is calculated as:")
print(f"c2.H = {S2} * {d} / {P} = {c2_H.numerator}/{c2_H.denominator}")
print("")
print("Finally, the Crawley-Nordström invariant c = (1/2) * c2.H is:")
print(f"c = (1/2) * ({c2_H.numerator}/{c2_H.denominator}) = {CN_invariant.numerator}/{CN_invariant.denominator}")