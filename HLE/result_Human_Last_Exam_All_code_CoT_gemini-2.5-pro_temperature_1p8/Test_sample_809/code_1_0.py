import math

# This script calculates the expected value of a game based on a series of probabilistic systems.
# We will solve for each probability step-by-step.

print("Step 1: Determine the value of h from system f.")
# Let P_f be the probability of successfully traversing f.
# P_f = x*1 + (1-x)*P_f^2, which is a quadratic equation: (1-x)P_f^2 - P_f + x = 0.
# Using the quadratic formula, the solutions for P_f are P_f = 1 and P_f = x/(1-x).
# Traversal is "guaranteed" (P_f=1) when the second solution is not a valid probability (<1).
# x/(1-x) >= 1  => x >= 1-x  => 2x >= 1  => x >= 0.5.
# h is the smallest value of x where this condition is met.
h = 0.5
print(f"The probability of traversing f is guaranteed (P(f)=1) when x >= 0.5.")
print(f"The smallest such value for x is h = {h}")

print("\nStep 2: Determine the traversal probability of g, P(g).")
# For g, P(direct) = h = 0.5. Let P(chain) be the probability of taking the recursive path.
# We are given P(hole) = 2 * P(chain).
# Probabilities sum to 1: h + P(hole) + P(chain) = 1
# 0.5 + 2*P(chain) + P(chain) = 1  =>  3*P(chain) = 0.5  => P(chain) = 1/6
prob_chain_g = 1.0 / 6.0
prob_direct_g = h
print(f"The probability of the direct path in g is h = {prob_direct_g}.")
print(f"The probability of the recursive path (chain of 6 g's) is 1/6.")
# The traversal probability P(g) satisfies: P(g) = P(direct)*1 + P(chain)*P(g)^6
# P(g) = 0.5 + (1/6) * P(g)^6.
# We solve this numerically using fixed-point iteration.
p_g = 0.5  # Initial guess
for _ in range(20):  # 20 iterations is more than enough for convergence
    p_g_next = prob_direct_g + prob_chain_g * (p_g ** 6)
    p_g = p_g_next
print(f"P(g) is the solution to P(g) = {prob_direct_g} + ({prob_chain_g:.4f}) * P(g)^6.")
print(f"Solving numerically, the probability of traversing g is P(g) = {p_g:.6f}")

print("\nStep 3: Determine the traversal probability of k, P(k).")
# k is a chain of 4 instances of g.
# P(k) = P(g)^4
p_k = p_g ** 4
print(f"k is a chain of 4 g's, so P(k) = P(g)^4.")
print(f"P(k) = ({p_g:.6f})^4 = {p_k:.6f}")

print("\nStep 4: Calculate the expected value of the game.")
n_trials = 100
p_success = p_k
win_threshold = 6 # We win if successes are >= 6
# We lose if successes are < 6 (i.e., 0, 1, 2, 3, 4, 5).
# Let's calculate P(lose) = P(successes < 6).
print(f"The game involves {n_trials} trials of k, each with a success probability p = {p_k:.6f}.")
print(f"We lose if the number of successes is less than {win_threshold}.")
prob_lose = 0.0
for i in range(win_threshold):
    # Binomial probability: C(n, i) * p^i * (1-p)^(n-i)
    term = math.comb(n_trials, i) * (p_success ** i) * ((1 - p_success) ** (n_trials - i))
    prob_lose += term
print(f"The probability of losing (fewer than {win_threshold} successes) is P(lose) = {prob_lose:.6f}.")

# Expected Value E = ($1 * P(win)) + (-$1 * P(lose))
# E = (1 * (1 - P(lose))) + (-1 * P(lose)) = 1 - 2 * P(lose)
expected_value = 1 - 2 * prob_lose
print("\nThe expected value equation is: E = (+1$) * P(win) - ($1) * P(lose)")
print(f"E = (1 - P(lose)) - P(lose) = 1 - 2 * P(lose)")
print(f"E = 1 - 2 * {prob_lose:.6f} = {expected_value:.6f}")
print("\nFinal Result:")
print(f"The expected value of playing one game is ${expected_value:.2f}.")