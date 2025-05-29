import math

# Parameters
n = 4  # number of trials
k = 3  # number of successes
p = 1/6  # probability of success on a single trial

# Binomial coefficient
binom_coeff = math.comb(n, k)

# Probability calculation
probability = binom_coeff * (p**k) * ((1-p)**(n-k))

print(probability)