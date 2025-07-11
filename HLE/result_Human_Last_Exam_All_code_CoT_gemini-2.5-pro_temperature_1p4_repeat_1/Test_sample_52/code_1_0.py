import math
from fractions import Fraction

# Step 1: Define the initial probabilities
prob_H = Fraction(1, 3)
prob_T = Fraction(2, 3)

# Let A be the event "all heads" (HHH).
# Let B be the event "even number of tails" (0 or 2 tails).

# Step 2: Calculate P(B), the probability of the condition.
# The probability of 0 tails (HHH)
prob_0_tails = prob_H ** 3

# The probability of 2 tails (e.g., HTT). There are C(3,2)=3 such combinations.
num_combinations_2_tails = math.comb(3, 2)
prob_2_tails = num_combinations_2_tails * (prob_H ** 1) * (prob_T ** 2)

# Total probability of an even number of tails
prob_B = prob_0_tails + prob_2_tails

# Step 3: Calculate P(A and B).
# The event "A and B" means all heads AND an even number of tails.
# This simplifies to just event A (all heads), since 0 tails is an even number.
# P(A and B) is the same as P(A), which is P(0 tails).
prob_A_and_B = prob_0_tails

# Step 4: Calculate the final conditional probability P(A|B) = P(A and B) / P(B)
final_probability = prob_A_and_B / prob_B

# Print the final calculation as an equation
print("The probability of all heads given an even number of tails is P(A|B) = P(A and B) / P(B)")
print("P(A and B) = P(All Heads) = (1/3)^3 = 1/27")
print("P(B) = P(0 Tails) + P(2 Tails) = 1/27 + 12/27 = 13/27")
print("\nFinal equation with numbers:")
print(f"({prob_A_and_B.numerator}/{prob_A_and_B.denominator}) / ({prob_B.numerator}/{prob_B.denominator}) = {final_probability.numerator}/{final_probability.denominator}")
print(f"\nThe probability of them all turning heads is: {final_probability}")
