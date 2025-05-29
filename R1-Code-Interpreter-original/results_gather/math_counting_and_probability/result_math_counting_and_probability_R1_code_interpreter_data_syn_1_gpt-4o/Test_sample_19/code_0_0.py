from math import comb

# Number of dice
n = 4
# Number of even numbers desired
k = 2
# Probability of rolling an even number on a single die
p = 1/2

# Calculate the probability using the binomial formula
probability = comb(n, k) * (p**k) * ((1-p)**(n-k))
print(probability)