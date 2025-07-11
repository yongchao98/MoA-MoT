from fractions import Fraction

# Step 1: Define the base probabilities for a single coin toss.
p_h = Fraction(1, 3)
p_t = Fraction(2, 3)

# The problem is a conditional probability: P(All Heads | Even Tails).
# The formula is P(All Heads and Even Tails) / P(Even Tails).

# The event "All Heads and Even Tails" is the same as "All Heads" (HHH),
# because having 0 tails is having an even number of tails.
# Let's calculate its probability.
p_all_heads = p_h * p_h * p_h

# Now, let's calculate the probability of the condition: "Even Tails".
# This can be 0 tails (HHH) or 2 tails (HTT, THT, TTH).

# We already have the probability for 0 tails.
p_0_tails = p_all_heads

# Calculate the probability for 2 tails. There are 3 combinations (HTT, THT, TTH).
# The probability of any one of these specific combinations (e.g., HTT) is:
p_one_combination_2_tails = p_h * p_t * p_t

# The total probability of getting 2 tails is 3 times this value.
p_total_2_tails = 3 * p_one_combination_2_tails

# The total probability of the condition "Even Tails" is the sum.
p_even_tails = p_0_tails + p_total_2_tails

# Finally, calculate the conditional probability.
conditional_probability = p_all_heads / p_even_tails

print("This problem is a conditional probability.")
print("We want to find P(All Heads | Number of Tails is Even).")
print("The formula is P(All Heads AND Even Tails) / P(Number of Tails is Even).\n")

print("First, the probability of 'All Heads' (HHH), which has 0 tails (an even number):")
print(f"P(HHH) = {p_h} * {p_h} * {p_h} = {p_all_heads}\n")

print("Next, the probability of having an even number of tails (0 or 2).")
print(f"P(0 Tails) is P(HHH) = {p_0_tails}")
print(f"P(2 Tails) is P(HTT, THT, TTH) = 3 * ({p_h} * {p_t}^2) = {p_total_2_tails}")
print(f"P(Even Tails) = P(0 Tails) + P(2 Tails) = {p_0_tails} + {p_total_2_tails} = {p_even_tails}\n")

print("Now we can find the final answer:")
print("P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)")
print(f"Final Equation: {p_all_heads} / {p_even_tails} = {conditional_probability}")
print(f"\nThe probability of them all turning heads is {conditional_probability}.")
<<<1/13>>>