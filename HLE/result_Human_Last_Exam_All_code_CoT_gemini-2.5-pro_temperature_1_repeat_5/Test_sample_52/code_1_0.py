from fractions import Fraction

# Step 1: Define the probabilities for a single coin flip.
p_h = Fraction(1, 3)
p_t = Fraction(2, 3)

# Step 2: Define the events.
# Event A: All three coins are heads (HHH). This implies 0 tails.
# Event B: The number of tails is even (0 or 2 tails).

# Step 3: Calculate the probability of the intersection of A and B.
# The event "A and B" is "all heads AND an even number of tails".
# If all coins are heads, the number of tails is 0, which is an even number.
# Therefore, the event "A and B" is the same as event "A".
# P(A and B) = P(A) = P(HHH)
p_A_and_B = p_h ** 3

# Step 4: Calculate the total probability of event B.
# Event B can happen in two ways: 0 tails or 2 tails.
# Probability of 0 tails (HHH):
p_0_tails = p_h ** 3
# Probability of 2 tails (HTT, THT, TTH):
# There are 3 combinations for 2 tails and 1 head.
# The probability of one such combination (e.g., HTT) is p_h * p_t**2
p_2_tails = 3 * p_h * (p_t ** 2)
# Total probability of event B:
p_B = p_0_tails + p_2_tails

# Step 5: Calculate the conditional probability P(A|B) = P(A and B) / P(B).
result = p_A_and_B / p_B

# Step 6: Print the explanation and the final equation.
print("This is a conditional probability problem: P(All Heads | Even Tails).")
print("The formula is P(A|B) = P(A and B) / P(B)\n")
print(f"The probability of the event 'All Heads and Even Tails', P(A and B), is the probability of HHH, which is (1/3)^3 = {p_A_and_B}.")
print(f"The probability of 'Even Tails' (0 or 2 tails), P(B), is P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails} = {p_B}.")
print("\nThe final equation is:")
print(f"{p_A_and_B} / {p_B} = {result}")
