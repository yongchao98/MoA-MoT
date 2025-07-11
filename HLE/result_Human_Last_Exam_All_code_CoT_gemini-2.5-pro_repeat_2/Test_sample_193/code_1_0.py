import math

def solve_poker_strategy():
    """
    This function determines and prints the optimal poker strategy for the given scenario.
    """
    # Game parameters
    pot = 10
    stack = 1000
    
    # As determined by the GTO analysis, the Expected Value (EV) for our value hand (AA)
    # is maximized by using the largest possible bet size. The relationship is:
    # EV(AA) = Pot + (Pot * BetSize) / (Pot + BetSize)
    # With a larger BetSize (B), the term (10*B)/(10+B) gets closer to 10, maximizing the EV.
    # Therefore, the optimal bet size is the entire stack.
    
    optimal_bet_size = stack

    # To satisfy the prompt, we output the numbers in the final EV equation.
    # This equation justifies our choice of bet sizing.
    print("The optimal strategy is derived from maximizing the Expected Value (EV) of betting with AA.")
    print("The EV for a bet of size B is: EV(AA) = pot + (pot * B) / (pot + B)")
    print(f"For the optimal all-in bet of {optimal_bet_size}, the equation with our numbers is:")
    ev_numerator = pot * optimal_bet_size
    ev_denominator = pot + optimal_bet_size
    final_ev = pot + ev_numerator / ev_denominator
    print(f"EV(AA) = {pot} + ({pot} * {optimal_bet_size}) / ({pot} + {optimal_bet_size}) = {final_ev:.2f}")

    # Since the all-in bet is the most profitable action with AA, we take it 100% of the time.
    action = "BET"
    sizing = optimal_bet_size
    frequency = 100

    # Round sizing and frequency to the nearest even number as per instructions.
    # In this case, 1000 and 100 are already even.
    rounded_sizing = int(round(sizing / 2.0)) * 2
    rounded_frequency = int(round(frequency / 2.0)) * 2
    
    # Final formatted strategy string
    final_answer_string = f"{action} {rounded_sizing} {rounded_frequency}%"
    
    print("\nThe optimal strategy for holding AA is therefore:\n")
    print(f"<<<{final_answer_string}>>>")

solve_poker_strategy()