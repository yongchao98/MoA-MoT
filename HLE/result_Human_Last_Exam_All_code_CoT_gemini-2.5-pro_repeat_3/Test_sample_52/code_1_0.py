from fractions import Fraction

# Step 1: Define the base probabilities for a single coin toss.
P_H = Fraction(1, 3)
P_T = Fraction(2, 3)

# Step 2: Define the events.
# Event A: All three coins are heads (HHH). This is also the case for 0 tails.
# Event B: The number of tails is even (0 or 2 tails).

# Step 3: Calculate the probability of the intersection of A and B.
# The event "All Heads" (A) is a case where there are 0 tails, which is an even number.
# Therefore, the event "A and B" is the same as event A.
# P(A and B) = P(HHH) = P(H) * P(H) * P(H)
prob_A_and_B = P_H ** 3

# Step 4: Calculate the probability of the condition (Event B).
# Event B occurs if there are 0 tails OR 2 tails.
# Probability of 0 tails is P(HHH), which we already calculated.
prob_0_tails = prob_A_and_B

# Probability of 2 tails. There are 3 combinations for this: HTT, THT, TTH.
# Each combination has a probability of P(H) * P(T) * P(T).
# So, P(2 tails) = 3 * P(H) * P(T)^2
prob_2_tails = 3 * P_H * (P_T ** 2)

# The total probability of Event B is the sum of these probabilities.
prob_B = prob_0_tails + prob_2_tails

# Step 5: Calculate the final conditional probability P(A|B).
# P(A|B) = P(A and B) / P(B)
final_prob = prob_A_and_B / prob_B

# Print the final equation with the calculated numbers
print("The probability of all heads given an even number of tails is calculated as:")
print("P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)")
print(f"P(All Heads | Even Tails) = P(0 Tails) / (P(0 Tails) + P(2 Tails))")
print(f"P(All Heads | Even Tails) = ({prob_0_tails}) / ({prob_0_tails} + {prob_2_tails})")
print(f"P(All Heads | Even Tails) = ({prob_0_tails}) / ({prob_B})")
print(f"The final probability is: {final_prob}")
print(f"As a decimal, this is approximately: {float(final_prob):.4f}")

<<<1/13>>>