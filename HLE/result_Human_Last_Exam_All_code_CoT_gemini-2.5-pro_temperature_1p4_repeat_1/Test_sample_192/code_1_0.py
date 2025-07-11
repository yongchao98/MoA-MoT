import math

def round_to_even(n):
    """Rounds a number to the nearest even integer."""
    return round(n / 2.0) * 2

def solve_poker_problem():
    """
    Solves the given poker problem by defining an optimal strategy
    and calculating the villain's response.
    """
    # Define an arbitrary optimal strategy. Any strategy where AA and QQ are played
    # identically is optimal, with an EV of $5. We choose one that requires rounding
    # to demonstrate the concept.
    
    # Strategy Choice
    bet_sizing_choice = 45
    bet_frequency_choice = 0.75  # 75%

    # Perform rounding as per instructions
    final_bet_sizing = round_to_even(bet_sizing_choice)
    final_bet_freq_pct = round_to_even(bet_frequency_choice * 100)
    final_check_freq_pct = 100 - final_bet_freq_pct

    # Since we must play AA and QQ identically to maintain an EV of $5,
    # the strategy is the same for both hands.
    action_string_A = f"BET ${final_bet_sizing} {final_bet_freq_pct}%, CHECK {final_check_freq_pct}%"
    action_string_B = action_string_A

    # Determine Villain's Response
    # As calculated in the analysis, villain's EV to call any bet is +$5,
    # while folding is $0. Thus, villain will always call.
    villain_response = "CALL 100%"

    # Format the final answer
    print(f"A) {action_string_A}")
    print(f"B) {action_string_B}")
    print(f"C) {villain_response}")

solve_poker_problem()