from fractions import Fraction

def calculate_strategy_probability(target_value, num_dice_to_reroll):
    """
    Calculates the probability of getting a Yahtzee for a given strategy.
    
    Args:
        target_value (int): The number we want all dice to show.
        num_dice_to_reroll (int): The number of dice we are rerolling.
        
    Returns:
        A Fraction object representing the probability of success.
    """
    # Probability of rolling the target value in a single throw
    p = Fraction(1, 2**target_value)
    
    # Probability of NOT rolling the target value
    q = 1 - p
    
    # With two rerolls, the probability of a single die NOT becoming the target value is q*q
    prob_single_die_fails = q**2
    
    # The probability of a single die successfully becoming the target value in up to two rolls
    prob_single_die_succeeds = 1 - prob_single_die_fails
    
    # The total probability is the success probability for one die raised to the power
    # of the number of dice we need to change, as they are independent events.
    total_prob = prob_single_die_succeeds**num_dice_to_reroll
    
    return total_prob

# --- Strategy A: Keep the '1' ---
# We keep one '1' and reroll the other 4 dice to try and get all '1's.
prob_A = calculate_strategy_probability(target_value=1, num_dice_to_reroll=4)

# --- Strategy B: Keep the three '3's ---
# We keep three '3's and reroll the other 2 dice to try and get all '3's.
prob_B = calculate_strategy_probability(target_value=3, num_dice_to_reroll=2)

# --- Calculate the Difference ---
# To subtract the fractions, we need a common denominator.
common_denominator = prob_A.denominator * prob_B.denominator // prob_A.denominator # A bit redundant but safe
prob_A_common = Fraction(prob_A.numerator * (prob_B.denominator // prob_A.denominator), prob_B.denominator)

# Calculate the difference
difference = prob_A - prob_B

print("Strategy A: Keep the '1' and aim for a Yahtzee of 1s.")
print(f"The probability of success is (3/4)^4 = {prob_A.numerator}/{prob_A.denominator}\n")

print("Strategy B: Keep the three '3's and aim for a Yahtzee of 3s.")
print(f"The probability of success is (15/64)^2 = {prob_B.numerator}/{prob_B.denominator}\n")

print("To find the difference, we use a common denominator:")
print(f"Probability A = {prob_A_common.numerator}/{prob_A_common.denominator}")
print(f"Probability B = {prob_B.numerator}/{prob_B.denominator}\n")

print("The difference in probabilities is:")
print(f"{prob_A_common.numerator}/{prob_A_common.denominator} - {prob_B.numerator}/{prob_B.denominator} = {difference.numerator}/{difference.denominator}")

# Final Answer
print(f"\n<<<{difference.numerator}/{difference.denominator}>>>")