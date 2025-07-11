import math

def solve_poker_strategy():
    """
    Calculates the optimal poker strategy based on the user's scenario.
    """
    pot = 10
    stack = 1000

    # The optimal bet sizing to exert maximum pressure with a polarized range is all-in.
    # The instruction is to round the sizing to the nearest even number. 1000 is already even.
    bet_size = 1000

    # Calculate the GTO bluffing frequency that makes the villain indifferent.
    # p_bluff_frequency = Bet / (Pot + Bet)
    p_bluff_frequency = bet_size / (pot + bet_size)
    p_bluff_percent = p_bluff_frequency * 100

    # Round the percentage to the nearest even number as per instructions.
    # The nearest even integer to 99.0099... is 100.
    # |99.0099 - 98| = 1.0099
    # |99.0099 - 100| = 0.9901
    # 100 is closer.
    rounded_bluff_percent = 100

    # Determine actions based on the calculated rounded frequency.
    aa_actions = f"BET ${bet_size} 100%"
    
    if rounded_bluff_percent == 100:
        # If we bluff 100% of the time with QQ, our strategy is identical for both hands.
        qq_actions = f"BET ${bet_size} 100%"
        # Villain's response to a bet from a 50/50 value/bluff range is to call.
        # EV(Call) = 0.5 * (Pot + Bet) + 0.5 * (-Bet) = 0.5 * Pot = $5.
        # EV(Fold) = $0. Villain calls.
        villain_response = "CALL 100%"
    else:
        # This case would handle other rounding results, leading to a different equilibrium.
        check_freq = 100 - rounded_bluff_percent
        qq_actions = f"BET ${bet_size} {rounded_bluff_percent}%, CHECK {check_freq}%"
        # In the equilibrium where we don't bluff 100%, villain is made indifferent and
        # chooses to fold to minimize our EV.
        villain_response = "FOLD 100%"

    # Print the final answer in the specified format.
    print(f"A) {aa_actions} B) {qq_actions} C) {villain_response}")

solve_poker_strategy()