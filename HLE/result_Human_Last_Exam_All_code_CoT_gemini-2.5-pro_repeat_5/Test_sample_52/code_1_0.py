import fractions

# This script calculates the conditional probability of getting all heads
# given that the number of tails from three coin flips is even.

# Step 1: Define the base probabilities for a single coin flip.
# Using the fractions module for precision.
p_heads = fractions.Fraction(1, 3)
p_tails = 1 - p_heads  # This will be 2/3

# Step 2: Define the events for the conditional probability formula P(A|B) = P(A and B) / P(B).
# Event A: All three coins are heads (HHH).
# Event B: The number of tails is even (0 or 2 tails).

# Step 3: Calculate P(A and B).
# This is the probability of getting all heads AND an even number of tails.
# Since "all heads" (HHH) has 0 tails, and 0 is an even number, this event
# is simply the probability of getting all heads.
p_A_and_B = p_heads * p_heads * p_heads

# Step 4: Calculate P(B), the probability of the condition being met.
# The number of tails can be even in two ways: 0 tails or 2 tails.
# Probability of 0 tails is the same as all heads.
p_0_tails = p_A_and_B
# Probability of 2 tails means 1 head and 2 tails. There are 3 such combinations (HTT, THT, TTH).
# The probability of one such combination (e.g., HTT) is p_heads * p_tails * p_tails.
p_2_tails = 3 * p_heads * (p_tails**2)
# The total probability of an even number of tails is the sum.
p_B = p_0_tails + p_2_tails

# Step 5: Calculate the final conditional probability, P(A|B).
final_probability = p_A_and_B / p_B

# Print the results, showing each number in the final equation.
print("This problem is solved using the conditional probability formula: P(A|B) = P(A and B) / P(B)\n")
print(f"The probability of our target event (All Heads), which satisfies the condition, is P(A and B) = {p_A_and_B}.")
print(f"The total probability of the condition (Even Tails) is P(B) = P(0 Tails) + P(2 Tails) = {p_0_tails} + {p_2_tails} = {p_B}.")
print("\nThe final equation is:")
print(f"P(All Heads | Even Tails) = P(A and B) / P(B) = {p_A_and_B} / {p_B} = {final_probability}")
print(f"\nThe final probability is {float(final_probability):.4f}")
