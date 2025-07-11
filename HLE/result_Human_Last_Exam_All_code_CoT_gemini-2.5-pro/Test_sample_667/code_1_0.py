from fractions import Fraction

# Step 1: Analyze Strategy A (Keep the '1')
# The goal is a Yahtzee of '1's. We need 4 more '1's.

# Probability of rolling a '1' on a single die
p1 = Fraction(1, 2)

# Probability for one die to become '1' within two re-rolls.
# This is P(success on 1st try) + P(fail on 1st try) * P(success on 2nd try)
q1 = p1 + (1 - p1) * p1

# Total probability for Strategy A is q1^4, as we need 4 dice to become '1'.
prob_A = q1**4

# Step 2: Analyze Strategy B (Keep the three '3's)
# The goal is a Yahtzee of '3's. We need 2 more '3's.

# Probability of rolling a '3' on a single die
p3 = Fraction(1, 8)

# Probability for one die to become '3' within two re-rolls.
q3 = p3 + (1 - p3) * p3

# Total probability for Strategy B is q3^2, as we need 2 dice to become '3'.
prob_B = q3**2

# Step 3: Calculate the difference
difference = prob_A - prob_B

# For a clear final equation, find a common denominator
common_denominator = prob_A.denominator * prob_B.denominator // difference.denominator
prob_A_common_num = prob_A.numerator * (common_denominator // prob_A.denominator)
prob_B_common_num = prob_B.numerator * (common_denominator // prob_B.denominator)

print("This script calculates the difference in probabilities between two Yahtzee strategies.")
print(f"Strategy A: Keep the '1'. The probability of success is P(A) = ({q1.numerator}/{q1.denominator})^4 = {prob_A.numerator}/{prob_A.denominator}")
print(f"Strategy B: Keep the three '3's. The probability of success is P(B) = ({q3.numerator}/{q3.denominator})^2 = {prob_B.numerator}/{prob_B.denominator}")
print("\nThe difference is P(A) - P(B):")
print(f"{prob_A.numerator}/{prob_A.denominator} - {prob_B.numerator}/{prob_B.denominator} = {prob_A_common_num}/{common_denominator} - {prob_B_common_num}/{common_denominator} = {difference.numerator}/{difference.denominator}")

# Final answer in the required format
final_answer = f"{difference.numerator}/{difference.denominator}"
print(f"\n<<< {final_answer} >>>")