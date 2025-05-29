import math

# Parameters
n = 5
k = 4
p = 2/3
q = 1 - p

# Binomial coefficient
binom_coeff = math.comb(n, k)

# Probability calculation
probability = binom_coeff * (p ** k) * (q ** (n - k))

# Output the result
print(probability)