from fractions import Fraction

# --- Introduction ---
print("This problem compares the probabilities of achieving a Yahtzee from the roll {1, 2, 3, 3, 3} using two different strategies.")
print("The probability of a single die showing 'k' is p(k) = 2**-k.")
print("We have two re-rolls available for any dice we don't keep.")
print("\nFirst, let's find the general probability of converting a single die to a target 'k' within two rolls.")
print("Let p = p(k). The probability of success is P_one_die = p (success on 1st try) + (1-p)*p (fail then succeed) = p*(2-p).")

# --- Strategy 1: Keep the '1' ---
print("\n--- Strategy 1: Keep the '1' ---")
print("We aim for five '1's, so we must convert the other 4 dice to '1's.")
p_1 = Fraction(1, 2)  # Probability of rolling a 1, p(1) = 2**-1
num_dice_1 = 4
prob_one_die_for_1 = p_1 * (2 - p_1)
prob_strategy_1 = prob_one_die_for_1 ** num_dice_1
print(f"The probability of turning one die into a '1' within two rolls is {p_1}*(2-{p_1}) = {prob_one_die_for_1}.")
print(f"The total probability of success for Strategy 1 is ({prob_one_die_for_1})^{num_dice_1}, which is {prob_strategy_1.numerator}/{prob_strategy_1.denominator}.")

# --- Strategy 2: Keep the three '3's ---
print("\n--- Strategy 2: Keep the three '3's ---")
print("We aim for five '3's, so we must convert the other 2 dice to '3's.")
p_3 = Fraction(1, 8)  # Probability of rolling a 3, p(3) = 2**-3
num_dice_2 = 2
prob_one_die_for_3 = p_3 * (2 - p_3)
prob_strategy_2 = prob_one_die_for_3 ** num_dice_2
print(f"The probability of turning one die into a '3' within two rolls is {p_3}*(2-{p_3}) = {prob_one_die_for_3}.")
print(f"The total probability of success for Strategy 2 is ({prob_one_die_for_3})^{num_dice_2}, which is {prob_strategy_2.numerator}/{prob_strategy_2.denominator}.")

# --- Final Calculation ---
print("\n--- Difference in Probabilities ---")
difference = prob_strategy_1 - prob_strategy_2

# To display the equation clearly, find a common denominator
p1_num = prob_strategy_1.numerator
p1_den = prob_strategy_1.denominator
p2_num = prob_strategy_2.numerator
p2_den = prob_strategy_2.denominator

# The common denominator is p2_den (4096)
p1_common_denom_num = p1_num * (p2_den // p1_den)

print(f"The probability of Strategy 1 is {p1_num}/{p1_den}.")
print(f"The probability of Strategy 2 is {p2_num}/{p2_den}.")
print(f"The difference is calculated as: ({p1_common_denom_num}/{p2_den}) - ({p2_num}/{p2_den}) = {difference.numerator}/{difference.denominator}")
print("\nFinal Answer:")
print(f"The numerical value of the difference is {float(difference)}.")
print("<<<1071/4096>>>")