import math

def solve_ur_probability():
    """
    Calculates and explains the probability of capturing a piece at square 12
    when starting with no pieces on the board in the Royal Game of Ur.
    """
    num_dice = 4
    total_outcomes = 2**num_dice
    target_square = 12
    rosette_square_on_entry = 4

    print("--- Royal Game of Ur Probability Analysis ---")
    print(f"Objective: Capture the opponent's piece on square {target_square}.")
    print("Your board state: No pieces currently on the board.")
    print("Rule: To enter the board, a roll of 'n' places a new piece on square 'n'.")
    print("-" * 45)

    # The number of ways to achieve each roll (0 to 4) is given by binomial coefficients C(4, k)
    # Roll 0: C(4,0) = 1 way
    # Roll 1: C(4,1) = 4 ways
    # Roll 2: C(4,2) = 6 ways
    # Roll 3: C(4,3) = 4 ways
    # Roll 4: C(4,4) = 1 way
    
    successful_outcomes = 0

    print("Analyzing possible moves for one turn:")
    # Scenario 1: Direct entry onto the target square.
    print(f"\n1. Can we land on square {target_square} directly?")
    print(f"   To do this, you must roll a {target_square}.")
    print(f"   The maximum possible roll with {num_dice} dice is {num_dice}.")
    print(f"   Result: Impossible. A roll of {target_square} is not possible.")

    # Scenario 2: Using the rosette rule.
    print(f"\n2. Can we use a rosette to reach square {target_square}?")
    print(f"   The only rosette you can land on when entering is square {rosette_square_on_entry}.")
    print(f"   This requires an initial roll of {rosette_square_on_entry}.")
    print(f"   If you roll a {rosette_square_on_entry}, you get a second roll.")
    needed_second_roll = target_square - rosette_square_on_entry
    print(f"   From square {rosette_square_on_entry}, you would need a second roll of {needed_second_roll} to reach {target_square}.")
    print(f"   The maximum possible roll is still {num_dice}.")
    print(f"   Result: Impossible. A second roll of {needed_second_roll} is not possible.")

    print("\n--- Conclusion ---")
    print("There are no possible dice rolls that allow you to move a piece to square 12 on your first turn.")
    
    # Final calculation
    numerator = successful_outcomes
    denominator = total_outcomes
    
    if numerator == 0:
        simplified_numerator = 0
        simplified_denominator = 1
    else:
        # This part is technically not needed since we know the numerator is 0,
        # but it's good practice for a general solution.
        common_divisor = math.gcd(numerator, denominator)
        simplified_numerator = numerator // common_divisor
        simplified_denominator = denominator // common_divisor

    print("\n--- Final Probability ---")
    print(f"Number of successful outcomes: {numerator}")
    print(f"Total number of possible outcomes: {denominator}")
    print(f"The final probability equation is: {numerator} / {denominator}")
    print(f"As a simplified common fraction, the probability is: {simplified_numerator}/{simplified_denominator}")

solve_ur_probability()