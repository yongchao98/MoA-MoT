import math

# Calculate binomial coefficients
def binomial_coefficient(n, k):
    return math.comb(n, k)

# Calculate the required values
binom_21_11 = binomial_coefficient(21, 11)
binom_20_10 = binomial_coefficient(20, 10)

# Find all n such that binom_20_n + binom_20_10 = binom_21_11
valid_n = []
for n in range(21):  # n can be from 0 to 20
    if binomial_coefficient(20, n) + binom_20_10 == binom_21_11:
        valid_n.append(n)

# Calculate the sum of all valid n
sum_valid_n = sum(valid_n)
print(sum_valid_n)