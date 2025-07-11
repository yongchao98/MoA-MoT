import math

# This script calculates the expected value of a game based on a series of probabilistic events.

# Step 1: Define h, the probability for the first path in g.
# Based on the analysis of structure f, the smallest value of x for which success is guaranteed is 0.5.
h = 0.5

# Step 2: Calculate the probability of successfully traversing g (p_g).
# The probability p_g is the solution to the equation: p_g = 0.5 + (1/6) * p_g^6.
# We solve this numerically using iteration.
p_g = 0.5  # Start with an initial guess.
# Iterate to find a stable solution for p_g.
for _ in range(10):
    p_g = 0.5 + (1/6) * (p_g ** 6)

# Step 3: Calculate the probability of successfully traversing k (p_k).
# k is a chain of four instances of g.
p_k = p_g ** 4

# Step 4: Calculate the expected value of the game.
# The number of successes S in n=100 trials follows a Binomial distribution.
# We use the Poisson approximation as n is large and p_k is small.
n = 100
p = p_k
lambda_val = n * p

# The opponent wins if S < 6 (i.e., S is 0, 1, 2, 3, 4, or 5).
# We calculate the opponent's win probability, P(S < 6).
prob_opponent_wins = 0
for i in range(6):  # This loops for i = 0, 1, 2, 3, 4, 5.
    poisson_term = (math.exp(-lambda_val) * (lambda_val ** i)) / math.factorial(i)
    prob_opponent_wins += poisson_term

# The expected value is given by the formula: EV = 1 - 2 * P(opponent wins).
# Here, we output the numbers that make up the final calculation as requested.
val1 = 1
val2 = 2
val3 = prob_opponent_wins
expected_value = val1 - val2 * val3

print(f"The probability of the opponent winning is P(S < 6) = {val3:.6f}")
print("The final equation for your expected value is: EV = A - B * C")
print(f"Where A = {val1}, B = {val2}, and C = P(S < 6) = {val3:.6f}")
print(f"EV = {val1} - {val2} * {val3:.6f} = {expected_value:.6f}")
print(f"Your expected value of playing one game, rounded to the nearest cent, is: ${expected_value:.2f}")
