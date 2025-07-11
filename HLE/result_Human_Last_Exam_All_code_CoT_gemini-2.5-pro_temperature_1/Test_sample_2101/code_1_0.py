import math

# This script calculates the difference in escape probabilities for an electron
# in an isosceles right triangle, based on the geometric derivation explained above.

# 1. Define the constant sqrt(2)
sqrt2 = math.sqrt(2)

# 2. Calculate the probability of escaping through the hypotenuse (P_H).
# This is derived as sqrt(2) / (2 + sqrt(2)).
prob_hypotenuse = sqrt2 / (2 + sqrt2)

# 3. Calculate the probability of escaping through a single leg (P_L).
# This is derived as 1 / (2 + sqrt(2)).
prob_one_leg = 1 / (2 + sqrt2)

# 4. The problem asks for the probability of escaping through "either of the two legs",
# which is the sum of the probabilities for each leg.
prob_both_legs = 2 * prob_one_leg

# 5. Calculate the final difference as requested.
difference = prob_hypotenuse - prob_both_legs

# 6. Output the components of the final calculation and the result.
# The final equation is: P(hypotenuse) - P(legs) = Difference
print("This script solves for the difference between electron escape probabilities.")
print(f"The probability of escaping through the hypotenuse is: {prob_hypotenuse:.5f}")
print(f"The probability of escaping through either of the two legs is: {prob_both_legs:.5f}")
print(f"The final equation is: {prob_hypotenuse:.5f} - {prob_both_legs:.5f} = {difference:.5f}")

# Also printing the exact symbolic result for clarity
# The difference simplifies to 2*sqrt(2) - 3
exact_difference = 2 * sqrt2 - 3
print(f"The exact value of the difference is 2*sqrt(2) - 3, which is approximately {exact_difference:.5f}")