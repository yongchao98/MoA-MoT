from fractions import Fraction

# Define the initial probabilities for Heads (H) and Tails (T)
p_H = Fraction(1, 3)
p_T = 1 - p_H

# Let A be the event of "all turning heads" (HHH).
# Let B be the event of "the number of tails is even" (0 or 2 tails).
# We want to find the conditional probability P(A|B) = P(A and B) / P(B).

print("--- Step-by-Step Calculation ---")

# Step 1: Calculate P(A and B), the probability of 'A and B' occurring together.
# The event 'A and B' is 'all heads' AND 'even number of tails'.
# The outcome 'HHH' has 0 tails, and 0 is an even number, so this satisfies both.
# P(A and B) is the probability of HHH.
p_A_and_B = p_H ** 3
print(f"1. The probability of getting all heads (which has 0 tails) is P(A and B) = (1/3)^3 = {p_A_and_B}")

# Step 2: Calculate P(B), the probability of the condition 'even number of tails'.
# This can be 0 tails (HHH) or 2 tails (HTT, THT, TTH).
# P(B) = P(0 tails) + P(2 tails)

# Probability of 0 tails is P(HHH)
p_0_tails = p_H ** 3

# Probability of 2 tails:
# There are 3 combinations for 2 tails (HTT, THT, TTH).
# The probability of one such combination (e.g., HTT) is P(H)*P(T)*P(T).
p_one_combo_2_tails = p_H * (p_T ** 2)
p_2_tails = 3 * p_one_combo_2_tails
p_B = p_0_tails + p_2_tails
print(f"2. The probability of an even number of tails is P(B) = P(0 tails) + P(2 tails) = {p_0_tails} + {p_2_tails} = {p_B}")

# Step 3: Calculate the final conditional probability P(A|B).
final_probability = p_A_and_B / p_B
print("\n--- Final Calculation ---")
print("The probability of them all turning heads GIVEN the number of tails is even is:")
print(f"P(All Heads | Even Tails) = P(A and B) / P(B)")
print(f"                            = {p_A_and_B} / {p_B}")
print(f"                            = {final_probability}")

<<<1/13>>>