import math

def solve_poker_strategy():
    """
    Calculates the optimal GTO strategy for the given river scenario.
    """
    # Define game parameters
    pot = 10.0

    # The theoretically optimal bet size (B) that maximizes our Expected Value (EV)
    # is B = pot / 2.
    optimal_b_unrounded = pot / 2

    # Per the user's request, we round the bet size to the nearest even number.
    # The nearest even numbers to 5 are 4 and 6. Both yield the same EV.
    # We will choose the larger bet size.
    bet_size = 6

    # A) Strategy with our value hand (AA):
    # We always bet for value.
    answer_A = f"BET ${int(bet_size)}"

    # B) Strategy with our bluffing hand (QQ):
    # We bluff with a frequency of B/P to make the villain indifferent.
    bluff_frequency = bet_size / pot
    bluff_percentage = int(bluff_frequency * 100)
    check_percentage = 100 - bluff_percentage

    answer_B = f"BET ${int(bet_size)} {bluff_percentage}%, CHECK {check_percentage}%"

    # C) Villain's response:
    # Our strategy makes the villain indifferent between calling and folding.
    # A "perfect nemesis" will choose the action that minimizes our EV.
    # Our EV is lower if the villain calls our bet than if they fold.
    # Therefore, the villain will always call.
    answer_C = "CALL"

    # Print the final formatted answer
    print(f"A) {answer_A} B) {answer_B} C) {answer_C}")

solve_poker_strategy()
<<<A) BET $6 B) BET $6 60%, CHECK 40% C) CALL>>>