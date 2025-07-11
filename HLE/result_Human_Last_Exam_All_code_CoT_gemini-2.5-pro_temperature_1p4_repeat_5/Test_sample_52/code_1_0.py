from fractions import Fraction

# Step 1: Define the base probabilities
p_h = Fraction(1, 3)  # Probability of a single head
p_t = Fraction(2, 3)  # Probability of a single tail

# Step 2: Calculate P(A and B)
# Event A is all heads (HHH). Event B is an even number of tails.
# The intersection 'A and B' is the event HHH, which has 0 tails (an even number).
# P(A and B) = P(HHH)
p_A_and_B = p_h ** 3

# Step 3: Calculate P(B)
# Event B is 0 tails or 2 tails.
# Probability of 0 tails is P(HHH)
p_0_tails = p_h ** 3

# Probability of 2 tails. There are 3 ways this can happen: HTT, THT, TTH.
# Each of these outcomes has a probability of p_h * p_t * p_t.
# Number of combinations for 2 tails from 3 coins is 3.
p_2_tails = 3 * (p_h ** 1) * (p_t ** 2)

# P(B) is the sum of the probabilities of these mutually exclusive cases.
p_B = p_0_tails + p_2_tails

# Step 4: Calculate the conditional probability P(A|B)
result = p_A_and_B / p_B

# --- Output the results ---
print("Let A be the event of all heads (HHH).")
print("Let B be the event of an even number of tails (0 or 2 tails).")
print("We need to calculate P(A|B) = P(A and B) / P(B).\n")

print(f"The probability of the intersection 'A and B' (which is just HHH) is:")
print(f"P(A and B) = (1/3)^3 = {p_A_and_B}\n")

print(f"The probability of B (even tails) is P(0 tails) + P(2 tails).")
print(f"P(0 tails) = P(HHH) = {p_0_tails}")
print(f"P(2 tails) = 3 * P(H) * P(T)^2 = 3 * {p_h} * ({p_t})^2 = {p_2_tails}")
print(f"P(B) = {p_0_tails} + {p_2_tails} = {p_B}\n")

print("The final conditional probability is:")
print(f"P(A|B) = P(A and B) / P(B) = {p_A_and_B} / {p_B} = {result}")
