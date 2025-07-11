import numpy as np
from scipy.optimize import fsolve
from scipy.stats import binom

# Step 1: Define h
h = 0.5

# Step 2: Solve for p_g, the probability of successfully traversing g.
# The equation is p_g^6 - 6*p_g + 3 = 0.
# We need to find the root that lies between 0 and 1.
func = lambda p: p**6 - 6 * p + 3
# An initial guess of 0.5 is reasonable as f(0.5) is close to 0.
p_g = fsolve(func, 0.5)[0]

# Step 3: Calculate p_k, the probability of successfully traversing k.
# k is a chain of 4 instances of g.
p_k = p_g**4

# Step 4: Calculate the expected value of the bet.
# Number of trials (traversal attempts of k)
n = 100
# The opponent wins if the number of successes is less than 6 (i.e., 0, 1, 2, 3, 4, or 5).
k_success_threshold = 5

# Probability of the opponent winning is P(S <= 5)
prob_opponent_wins = binom.cdf(k_success_threshold, n, p_k)

# The final equation for expected value is EV = 1 - 2 * P(Opponent Wins)
expected_value = 1 - 2 * prob_opponent_wins

print(f"The value of h is: {h}")
print(f"The success probability of g (p_g) is: {p_g:.6f}")
print(f"The success probability of k (p_k) is: {p_k:.6f}")
print("\nCalculating the expected value:")
print(f"The opponent wins if successes are less than 6. The probability of this is P(S <= 5).")
print(f"P(S <= 5) = {prob_opponent_wins:.6f}")
print(f"The final equation is: EV = 1 - 2 * {prob_opponent_wins:.6f}")
print(f"Your expected value is: ${expected_value:.2f}")
