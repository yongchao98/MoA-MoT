import math
from fractions import Fraction

# Define the probabilities for a single coin toss
p_h = Fraction(1, 3)  # Probability of Heads
p_t = Fraction(2, 3)  # Probability of Tails
n_coins = 3

# Calculate P(A and B), where A is "all heads" and B is "even number of tails".
# Since "all heads" (0 tails) is an event with an even number of tails, P(A and B) is simply P(A).
# P(A) = P(HHH) = P(H)^3
prob_A_and_B = p_h ** n_coins

# Calculate P(B), the probability of an even number of tails (0 or 2 tails).
# P(0 tails) = P(HHH)
prob_zero_tails = p_h ** 3

# P(2 tails) = (Number of ways to get 2 tails) * P(T)^2 * P(H)^1
# The number of ways is C(3, 2) = 3.
num_comb_two_tails = math.comb(n_coins, 2)
prob_two_tails = num_comb_two_tails * (p_t ** 2) * (p_h ** 1)

# P(B) = P(0 tails) + P(2 tails)
prob_B = prob_zero_tails + prob_two_tails

# Calculate the final conditional probability P(A|B) = P(A and B) / P(B)
final_probability = prob_A_and_B / prob_B

# Print the components of the final equation as requested
print("The probability of event 'A and B' (all heads, which has an even number of tails) is P(A and B) = P(HHH).")
print(f"P(A and B) = ({p_h})^3 = {prob_A_and_B}")
print("\nThe probability of event 'B' (even number of tails) is P(0 tails) + P(2 tails).")
print(f"P(B) = {prob_zero_tails} + {prob_two_tails} = {prob_B}")
print("\nThe final conditional probability P(A|B) is calculated as P(A and B) / P(B).")
print("\nFinal Equation:")
print(f"{prob_A_and_B} / {prob_B} = {final_probability}")

<<<1/13>>>