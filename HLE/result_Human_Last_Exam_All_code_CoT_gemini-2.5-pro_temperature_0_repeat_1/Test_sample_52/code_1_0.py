from fractions import Fraction
import math

# Define the probabilities for a single coin flip
p_heads = Fraction(1, 3)
p_tails = Fraction(2, 3)
n_coins = 3

# The event of interest is "all heads", which means 0 tails.
# Using the binomial probability formula: C(n, k) * p^k * (1-p)^(n-k)
# Let's define getting "Tails" as a success.
# Probability of 0 tails (0 successes):
k_0_tails = 0
prob_0_tails = math.comb(n_coins, k_0_tails) * (p_tails**k_0_tails) * (p_heads**(n_coins - k_0_tails))

# The condition is an even number of tails (0 or 2).
# Probability of 2 tails (2 successes):
k_2_tails = 2
prob_2_tails = math.comb(n_coins, k_2_tails) * (p_tails**k_2_tails) * (p_heads**(n_coins - k_2_tails))

# Total probability of the condition being met (even number of tails)
prob_even_tails = prob_0_tails + prob_2_tails

# The event "all heads" is the same as "0 tails".
prob_all_heads = prob_0_tails

# The conditional probability is P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)
final_probability = prob_all_heads / prob_even_tails

# Print the final equation with each number
print("The probability of all heads (0 tails) is P(A) = {}".format(prob_all_heads))
print("The probability of an even number of tails is P(B) = P(0 tails) + P(2 tails) = {} + {} = {}".format(prob_0_tails, prob_2_tails, prob_even_tails))
print("\nThe final conditional probability is P(A|B) = P(A) / P(B)")
print("Final Equation: {} / {} = {}".format(prob_all_heads, prob_even_tails, final_probability))
<<<1/13>>>