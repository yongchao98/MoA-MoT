from fractions import Fraction

# Step 1: Define a function to calculate the probability of a single die
# becoming a target number within two rolls.
def get_prob_one_die(target_number):
    """
    Calculates the probability of a single die becoming the target number
    within two available rolls.
    The probability of rolling the target number is p = 2**(-target_number).
    The chance of success is p + (1-p)*p.
    """
    p = Fraction(1, 2**target_number)
    prob_success = p + (1 - p) * p
    return prob_success

# Step 2: Calculate the total probability for Strategy A (keeping the '1').
# We need to turn 4 dice into '1's.
target_A = 1
num_dice_A = 4
prob_one_die_A = get_prob_one_die(target_A)
prob_A = prob_one_die_A ** num_dice_A

# Step 3: Calculate the total probability for Strategy B (keeping the three '3's).
# We need to turn 2 dice into '3's.
target_B = 3
num_dice_B = 2
prob_one_die_B = get_prob_one_die(target_B)
prob_B = prob_one_die_B ** num_dice_B

# Step 4: Calculate the difference between the two probabilities.
difference = prob_A - prob_B

# Step 5: Print the components of the final equation and the result.
print("This script calculates the difference between two Yahtzee strategies.")
print("-" * 20)
print(f"Strategy A: Keep the '1' and aim for five '1's.")
print(f"The probability of success is P(A) = {prob_A.numerator}/{prob_A.denominator}")
print("-" * 20)
print(f"Strategy B: Keep the three '3's and aim for five '3's.")
print(f"The probability of success is P(B) = {prob_B.numerator}/{prob_B.denominator}")
print("-" * 20)
print("The difference between the probabilities is:")
print(f"P(A) - P(B) = {prob_A.numerator}/{prob_A.denominator} - {prob_B.numerator}/{prob_B.denominator} = {difference.numerator}/{difference.denominator}")
print(f"The numerical value of the difference is approximately {float(difference):.6f}")
