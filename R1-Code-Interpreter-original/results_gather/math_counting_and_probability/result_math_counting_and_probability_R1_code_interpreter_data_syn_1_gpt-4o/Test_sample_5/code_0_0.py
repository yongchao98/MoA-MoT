from math import comb

n = 9
k = 1
a = 1  # Coefficient of x
b = -1  # Coefficient of -1

# Calculate the binomial coefficient
binomial_coefficient = comb(n, k)

# Calculate the term
term_coefficient = binomial_coefficient * (a ** (n - k)) * (b ** k)

print(term_coefficient)