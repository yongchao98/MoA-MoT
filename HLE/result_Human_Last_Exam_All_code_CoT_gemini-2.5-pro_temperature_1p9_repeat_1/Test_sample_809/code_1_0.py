import numpy as np
from scipy.optimize import fsolve
from scipy.stats import binom

# Step 1: Find the value of h for structure f.
# The probability of traversing f, P(f), is given by the recursive formula:
# P(f) = x * 1 + (1 - x) * P(f)^2
# We can rearrange this into a quadratic equation: (1 - x) * P(f)^2 - P(f) + x = 0.
# We need to find the smallest value of x, which we call h, that "guarantees" traversal, meaning P(f) = 1.
# Let's substitute P(f)=1 into the quadratic equation to see when this holds:
# (1 - x) * (1)^2 - 1 + x = 1 - x - 1 + x = 0.
# This identity holds for any x, which means P(f)=1 is always a mathematical solution.
# However, the physically meaningful probability is the smaller of the two solutions of the quadratic formula.
# The solutions are P(f) = 1 and P(f) = x / (1 - x).
# For the probability to be guaranteed (P(f)=1), the other solution must not be a valid, smaller probability.
# The crossover point is when x / (1 - x) = 1, which solves to 2x = 1, or x = 0.5.
# For x >= 0.5, the only valid probability solution is P(f)=1.
# Therefore, the smallest value of x that guarantees traversal is 0.5.
h = 0.5

# Step 2: Determine probabilities for g.
# The three paths in g have probabilities p_direct, p_hole, and p_chain.
# p_direct = h = 0.5
# p_direct + p_hole + p_chain = 1
# p_hole = 2 * p_chain
# Substituting gives: 0.5 + 2 * p_chain + p_chain = 1 => 3 * p_chain = 0.5
p_chain = 0.5 / 3
p_hole = 2 * p_chain
p_direct = h

# Step 3: Solve for the traversal probability of g, P(g).
# P(g) follows a recursive formula: P(g) = p_direct * 1 + p_hole * 0 + p_chain * P(g)^6
def g_equation(p_g):
    # This is the equation: P(g) = 0.5 + (1/6) * P(g)^6
    return p_g - p_direct - p_chain * p_g**6

# We solve this equation numerically. The solution will be close to 0.5.
initial_guess = 0.5
p_g = fsolve(g_equation, initial_guess)[0]

# Step 4: Calculate the traversal probability of k, P(k).
# k is a chain of four g's.
p_k = p_g**4

# Step 5: Calculate the expected value of the game.
n_trials = 100
# The opponent wins if the number of successes, N, is less than 6 (i.e., N <= 5).
k_successes = 5
prob_lose_bet = binom.cdf(k_successes, n_trials, p_k)

# Expected Value (EV) = P(win) * ($1) + P(lose) * (-$1)
# EV = (1 - prob_lose_bet) - prob_lose_bet
expected_value = 1 - 2 * prob_lose_bet

print(f"The value of h is {h}")
print(f"The probabilities for the paths in g are:")
print(f"  - Direct path: {p_direct}")
print(f"  - Path to hole: {p_hole:.4f}")
print(f"  - Recursive chain path: {p_chain:.4f}")
print("")
print(f"The traversal probability of g, P(g), is the solution to P(g) = {p_direct} + {p_chain:.4f} * P(g)^6")
print(f"Solving this gives P(g) ~= {p_g:.6f}")
print("")
print(f"The traversal probability of k is P(k) = P(g)^4 = ({p_g:.6f})^4 ~= {p_k:.6f}")
print("")
print("For the game:")
print(f"We need the probability of less than 6 successes in {n_trials} trials, with success probability p = {p_k:.6f}")
print(f"The probability of you losing the bet, P(N < 6), is {prob_lose_bet:.6f}")
print("")
print(f"The expected value equation is: EV = $1 * (1 - P(N<6)) + (-$1) * P(N<6)")
print(f"Final calculation: EV = 1 - 2 * {prob_lose_bet:.4f} = {expected_value:.4f}")
print(f"\nYour expected value of playing one game, rounded to the nearest cent, is: ${expected_value:.2f}")
<<<O>>>