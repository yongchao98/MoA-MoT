from math import comb

# Number of dice
n = 6
# Number of even numbers needed
k = 3
# Probability of rolling an even number
p = 0.5

# Calculate the probability using the binomial formula
probability = comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

print(probability)