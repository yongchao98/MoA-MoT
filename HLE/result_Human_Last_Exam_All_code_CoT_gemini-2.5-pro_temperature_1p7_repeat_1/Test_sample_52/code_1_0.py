from fractions import Fraction

# Define the basic probabilities for a single coin toss
P_H = Fraction(1, 3)
P_T = Fraction(2, 3)

# Define Event A: All 3 coins are heads (HHH)
# Define Event B: The number of tails is even (0 or 2 tails)

# We want to find P(A | B) = P(A and B) / P(B)

# Calculate P(A and B), the probability of 'all heads' and 'even tails'.
# This simplifies to just P(all heads), since 0 tails is an even number.
prob_A_and_B = P_H**3

# Calculate P(B), the probability of an even number of tails.
# This is P(0 tails) + P(2 tails).
prob_0_tails = P_H**3
# For 2 tails, there are 3 combinations (HTT, THT, TTH).
# Each has probability P_H * P_T**2.
prob_2_tails = 3 * (P_H * P_T**2)
prob_B = prob_0_tails + prob_2_tails

# Calculate the final conditional probability P(A | B).
final_prob = prob_A_and_B / prob_B

# Print the final equation with all the numbers.
print("The conditional probability of all coins turning heads, given an even number of tails, is calculated as:")
print("P(all heads | even tails) = P(all heads and even tails) / P(even tails)")
print(f"P(all heads | even tails) = P(0 tails) / (P(0 tails) + P(2 tails))")
print(f"P(all heads | even tails) = ({prob_A_and_B.numerator}/{prob_A_and_B.denominator}) / (({prob_0_tails.numerator}/{prob_0_tails.denominator}) + ({prob_2_tails.numerator}/{prob_2_tails.denominator}))")
print(f"P(all heads | even tails) = ({prob_A_and_B.numerator}/{prob_A_and_B.denominator}) / ({prob_B.numerator}/{prob_B.denominator})")
print(f"P(all heads | even tails) = {final_prob.numerator}/{final_prob.denominator}")
