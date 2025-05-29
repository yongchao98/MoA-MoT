import math

n = 9
k = 3
binomial_coefficient = math.comb(n, k)
constant_term = binomial_coefficient * (5 ** k)

print(constant_term)