from fractions import Fraction

# Step 1: Define the base probabilities
p_H = Fraction(1, 3)
p_T = 1 - p_H

# Step 2 & 4: Define the target event (A) and the conditional event (B)
# Event A: All heads (HHH). This is also the event of 0 tails.
p_A = p_H**3

# Event B: The number of tails is even (0 or 2 tails).
# Probability of 0 tails is the same as P(A).
p_0_tails = p_A
# Probability of 2 tails. There are 3 combinations (HTT, THT, TTH),
# each with probability P(H) * P(T)^2.
p_2_tails = 3 * (p_H**1) * (p_T**2)
# Total probability of the condition (Event B).
p_B = p_0_tails + p_2_tails

# Step 5: Calculate the conditional probability P(A|B)
# The event "A and B" is "All Heads AND Even Tails". Since "All Heads" (0 tails)
# is a case of "Even Tails", the intersection is just "All Heads" (A).
# So, P(A and B) = P(A).
p_A_given_B = p_A / p_B

# --- Output the results ---
print("This script calculates the probability of getting all heads given that the number of tails is even.")
print(f"\nProbability of a single Head P(H) = {p_H}")
print(f"Probability of a single Tail P(T) = {p_T}")

print("\nLet A be the event 'all heads' (HHH).")
print(f"P(A) = P(H)^3 = ({p_H})^3 = {p_A}")

print("\nLet B be the event 'the number of tails is even' (0 or 2 tails).")
print(f"P(0 Tails) = P(HHH) = {p_0_tails}")
print(f"P(2 Tails) = 3 * P(H) * P(T)^2 = 3 * {p_H} * ({p_T})^2 = {p_2_tails}")
print(f"P(B) = P(0 Tails) + P(2 Tails) = {p_0_tails} + {p_2_tails} = {p_B}")

print("\nWe want to find P(A | B) = P(A and B) / P(B).")
print("Since event A is a subset of event B, P(A and B) = P(A).")
print("\nThe final equation is:")
print(f"P(A | B) = P(A) / P(B) = ({p_A}) / ({p_B}) = {p_A_given_B}")
print(f"In numbers: ({p_A.numerator}/{p_A.denominator}) / ({p_B.numerator}/{p_B.denominator}) = {p_A_given_B.numerator}/{p_A_given_B.denominator}")
