import fractions

# Use the fractions module for precise fractional arithmetic
p_h = fractions.Fraction(1, 3)
p_t = fractions.Fraction(2, 3)

# 1. Calculate the probability of the event we are interested in: all heads (HHH).
# This is also the probability of having 0 tails.
p_all_heads = p_h ** 3

# 2. Calculate the probability of the condition: an even number of tails (0 or 2).
# The probability of 0 tails is p_all_heads, which we already calculated.
p_0_tails = p_all_heads

# The probability of 2 tails. There are 3 combinations for this: TTH, THT, HTT.
# Each combination has the same probability: p_t * p_t * p_h
p_one_combination_2_tails = (p_t ** 2) * (p_h ** 1)
p_2_tails = 3 * p_one_combination_2_tails

# The total probability of the condition (even number of tails) is P(0 tails) + P(2 tails).
p_condition_even_tails = p_0_tails + p_2_tails

# 3. Calculate the conditional probability.
# P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)
# Note: P(All Heads and Even Tails) is the same as P(All Heads) because
# the 'All Heads' event is a subset of the 'Even Tails' event.
conditional_prob = p_all_heads / p_condition_even_tails

# 4. Print the results and the final equation.
print("Step 1: Calculate the probability of 'all heads' (which implies 0 tails).")
print(f"P(All Heads) = (1/3)^3 = {p_all_heads}\n")

print("Step 2: Calculate the probability of the condition 'even number of tails'.")
print("This means 0 tails or 2 tails.")
print(f"P(0 Tails) = P(All Heads) = {p_0_tails}")
print(f"P(2 Tails) = 3 * (1/3)^1 * (2/3)^2 = {p_2_tails}")
print(f"P(Even Tails) = P(0 Tails) + P(2 Tails) = {p_0_tails} + {p_2_tails} = {p_condition_even_tails}\n")

print("Step 3: Calculate the conditional probability.")
print("P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)\n")

print("Final equation with the calculated probabilities:")
# The following line prints each number in the final equation.
print(f"{p_all_heads} / {p_condition_even_tails} = {conditional_prob}")
