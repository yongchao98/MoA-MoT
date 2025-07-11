import math
from fractions import Fraction

def calculate_yahtzee_probability(target_k, num_kept):
    """
    Calculates the total probability of achieving a Yahtzee of 'target_k'
    over two re-rolls, starting with 'num_kept' dice already showing 'target_k'.

    Args:
        target_k (int): The number we are aiming for (e.g., 1 or 3).
        num_kept (int): The number of dice we are keeping that show the target_k.

    Returns:
        Fraction: The total probability of success.
    """
    N = 5  # Total number of dice

    # Probability of rolling the target number k
    p = Fraction(1, 2**target_k)
    # Probability of NOT rolling the target number k
    q = 1 - p

    # --- First re-roll ---
    num_to_roll_1 = N - num_kept
    # Probability of success on the first re-roll (getting all remaining dice to be target_k)
    prob_succ_roll_1 = p**num_to_roll_1

    # --- Second re-roll ---
    # This is the probability of succeeding on the 2nd roll, GIVEN failure on the 1st.
    prob_succ_on_roll_2_after_failure = Fraction(0)

    # Iterate through the possible number of new 'k's (j) on the first re-roll.
    # Since we assume failure on the first roll, j is less than the number of dice rolled.
    for j in range(num_to_roll_1):
        # Probability of getting exactly j new 'k's from the first re-roll
        # This uses the binomial probability formula: C(n, k) * p^k * (1-p)^(n-k)
        prob_of_getting_j_new_ks = (
            math.comb(num_to_roll_1, j) *
            (p**j) *
            (q**(num_to_roll_1 - j))
        )

        # If we got j new k's, we now have this many k's in total for the second roll.
        num_k_for_roll_2 = num_kept + j
        num_to_roll_2 = N - num_k_for_roll_2

        # Probability of succeeding on the second re-roll given this state
        prob_succ_given_j_new_ks = p**num_to_roll_2

        # Add the contribution of this specific path to the total
        prob_succ_on_roll_2_after_failure += prob_of_getting_j_new_ks * prob_succ_given_j_new_ks

    # The total probability is succeeding on roll 1 OR failing roll 1 and succeeding on roll 2.
    total_prob = prob_succ_roll_1 + prob_succ_on_roll_2_after_failure
    return total_prob

# Calculate probability for Strategy 1: Keep one '1'
prob_keep_1 = calculate_yahtzee_probability(target_k=1, num_kept=1)

# Calculate probability for Strategy 2: Keep three '3's
prob_keep_3s = calculate_yahtzee_probability(target_k=3, num_kept=3)

# Find the common denominator to display the fractions nicely
common_denominator = math.lcm(prob_keep_1.denominator, prob_keep_3s.denominator)

# Numerators for the common denominator
num1 = prob_keep_1.numerator * (common_denominator // prob_keep_1.denominator)
num2 = prob_keep_3s.numerator * (common_denominator // prob_keep_3s.denominator)

# Calculate the difference
difference = prob_keep_1 - prob_keep_3s
diff_num = difference.numerator
diff_den = difference.denominator

print("Chance of success if you keep the '1':")
print(f"P(keep 1) = {num1}/{common_denominator}")

print("\nChance of success if you keep the three '3's:")
print(f"P(keep 3s) = {num2}/{common_denominator}")

print("\nThe difference between the two chances is:")
print(f"{num1}/{common_denominator} - {num2}/{common_denominator} = {diff_num}/{diff_den}")

print(f"\n<<<1071/4096>>>")