from fractions import Fraction

# This script calculates the difference in probabilities for two Yahtzee strategies.

# --- Strategy A: Keep the '1' and aim for five 1s ---

# The target number is 1. We need to change 4 dice.
target_a = 1
dice_to_change_a = 4

# The probability of rolling a 1 is p(1) = 2^(-1).
p_a = Fraction(1, 2**target_a)

# The probability of a single die becoming a 1 within two rolls is:
# P(success on 1st roll) + P(fail on 1st roll) * P(success on 2nd roll)
prob_one_die_a = p_a + (1 - p_a) * p_a

# The total probability for Strategy A is the single-die probability to the power of the number of dice to change.
prob_total_a = prob_one_die_a ** dice_to_change_a

# --- Strategy B: Keep the three '3's and aim for five 3s ---

# The target number is 3. We need to change 2 dice.
target_b = 3
dice_to_change_b = 2

# The probability of rolling a 3 is p(3) = 2^(-3).
p_b = Fraction(1, 2**target_b)

# The probability of a single die becoming a 3 within two rolls.
prob_one_die_b = p_b + (1 - p_b) * p_b

# The total probability for Strategy B.
prob_total_b = prob_one_die_b ** dice_to_change_b

# --- Calculate and Print the Difference ---

# The difference between the probabilities of the two strategies.
difference = prob_total_a - prob_total_b

print("This program calculates the difference between two Yahtzee strategies.")

print("\n--- Strategy A: Aiming for five 1s ---")
print(f"The probability of a die becoming a 1 in two rolls is ({p_a}) + (1 - {p_a}) * ({p_a}) = {prob_one_die_a}")
print(f"The total probability P(A) is ({prob_one_die_a})^4 = {prob_total_a}")

print("\n--- Strategy B: Aiming for five 3s ---")
print(f"The probability of a die becoming a 3 in two rolls is ({p_b}) + (1 - {p_b}) * ({p_b}) = {prob_one_die_b}")
print(f"The total probability P(B) is ({prob_one_die_b})^2 = {prob_total_b}")

# To show the calculation with a common denominator
prob_a_common_denom = prob_total_a.limit_denominator(prob_total_b.denominator)

print("\n--- Final Difference ---")
print(f"The difference is P(A) - P(B)")
print(f"= {prob_total_a} - {prob_total_b}")
print(f"= {prob_a_common_denom} - {prob_total_b}")
print(f"= {difference.numerator}/{difference.denominator}")
print(f"\nThe final numerical difference is: {float(difference)}")
