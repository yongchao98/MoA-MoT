from fractions import Fraction

# Step 1: Define the probability for a single coin toss.
# P(H) is the probability of heads, P(T) is the probability of tails.
p_h = Fraction(1, 3)
p_t = 1 - p_h

print(f"The probability of a single coin turning heads is P(H) = {p_h}")
print(f"The probability of a single coin turning tails is P(T) = {p_t}")
print("\nIt is known that the number of tails is even. This can happen in two ways:")
print("1. Zero tails (HHH)")
print("2. Two tails (HTT, THT, TTH)")

# Step 2: Calculate the probability of the outcome with zero tails (HHH).
p_hhh = p_h * p_h * p_h

print(f"\nThe probability of zero tails (HHH) is: {p_h} * {p_h} * {p_h} = {p_hhh}")

# Step 3: Calculate the probability of having two tails.
# There are 3 combinations for two tails (HTT, THT, TTH).
# Each has a probability of p_h * p_t * p_t.
p_2_tails = 3 * (p_h * p_t * p_t)

print(f"The probability of two tails is: 3 * ({p_h} * {p_t} * {p_t}) = {p_2_tails}")

# Step 4: Calculate the total probability of the condition (even number of tails).
# This is the sum of the probabilities from steps 2 and 3.
p_condition_met = p_hhh + p_2_tails

print(f"\nThe total probability of the number of tails being even is: {p_hhh} + {p_2_tails} = {p_condition_met}")

# Step 5: Calculate the conditional probability.
# This is P(HHH) / P(even tails).
final_probability = p_hhh / p_condition_met

print("\nThe probability of them all turning heads given the number of tails is even is:")
print(f"P(all heads | even tails) = P(all heads) / P(even tails)")
print(f"                          = {p_hhh} / {p_condition_met}")
print(f"                          = {final_probability}")

<<<1/13>>>