from fractions import Fraction

def calculate_yahtzee_prob(kept_dice_count, kept_dice_value):
    """
    Calculates the probability of getting a Yahtzee given the number of dice kept
    and their value, with two rolls remaining.
    The formula used is (2p - p^2)^m, where:
    p = probability of rolling the desired number (2^-n)
    m = number of dice to re-roll (5 - k)
    """
    m = 5 - kept_dice_count
    p = Fraction(1, 2**kept_dice_value)
    prob = (2 * p - p**2)**m
    return prob

# --- Strategy A: Keep the single '1' ---
prob_A = calculate_yahtzee_prob(kept_dice_count=1, kept_dice_value=1)

# --- Strategy B: Keep the three '3's ---
prob_B = calculate_yahtzee_prob(kept_dice_count=3, kept_dice_value=3)

# --- Calculate the difference ---
difference = prob_A - prob_B

# --- Print the detailed results ---
print("This script calculates the difference in Yahtzee probability between two strategies.")
print("\n--- Strategy A: Keep the single '1' ---")
print("We keep 1 die and need 4 more '1's.")
print("The probability of rolling a '1' is p = 1/2.")
print(f"The total probability of success P(A) is (3/4)^4 = {prob_A}")

print("\n--- Strategy B: Keep the three '3's ---")
print("We keep 3 dice and need 2 more '3's.")
print("The probability of rolling a '3' is p = 1/8.")
print(f"The total probability of success P(B) is (15/64)^2 = {prob_B}")

print("\n--- Difference Calculation ---")
# To show the common denominator calculation, we find the denominator of the larger fraction
common_denominator = prob_B.denominator
prob_A_common_denom = prob_A.limit_denominator(common_denominator)

print(f"Difference = P(A) - P(B)")
print(f"Difference = {prob_A} - {prob_B}")
print(f"Difference = {prob_A_common_denom} - {prob_B}")
print(f"Difference = {difference}")
print(f"\nAs a decimal, the difference is approximately {float(difference):.6f}")

# Final answer in the required format
final_answer = float(difference)
# print(f"\n<<<{final_answer}>>>")