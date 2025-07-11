from fractions import Fraction

# Step 1: Define the base probabilities for a single coin flip.
p_heads = Fraction(1, 3)
p_tails = Fraction(2, 3)

# Step 2: Define the event A: all three coins are heads (HHH).
# This is also the event of having 0 tails.
p_event_A = p_heads ** 3

# Step 3: Define the event B: the number of tails is even (0 or 2 tails).
# Calculate the probability of 0 tails.
p_0_tails = p_event_A
print(f"The probability of 0 tails (HHH) is (1/3)^3 = {p_0_tails}")

# Calculate the probability of 2 tails.
# There are 3 possible combinations for 2 tails and 1 head (HTT, THT, TTH).
# The probability for any single one of these is (1/3) * (2/3)^2.
num_combinations_2_tails = 3
p_2_tails = num_combinations_2_tails * p_heads * (p_tails ** 2)
print(f"The probability of 2 tails is 3 * (1/3) * (2/3)^2 = {p_2_tails}")

# The total probability of event B is the sum of these probabilities.
p_event_B = p_0_tails + p_2_tails
print(f"The total probability of an even number of tails is P(B) = {p_0_tails} + {p_2_tails} = {p_event_B}")
print("-" * 20)

# Step 4: Calculate the conditional probability P(A|B) = P(A and B) / P(B).
# Since event A (all heads) implies an even number of tails (0), the event "A and B" is just A.
# Therefore, P(A|B) = P(A) / P(B).
p_A_given_B = p_event_A / p_event_B

# Step 5: Print the final calculation and result.
print("The desired probability is P(All Heads | Even Tails).")
print("This is calculated as P(All Heads) / P(Even Tails).")
print("\nThe final equation is:")
print(f"{p_event_A} / {p_event_B} = {p_A_given_B}")
<<<1/13>>>