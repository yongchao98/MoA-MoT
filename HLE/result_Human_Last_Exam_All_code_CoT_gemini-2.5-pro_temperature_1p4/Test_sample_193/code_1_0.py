import math

def solve_poker_strategy():
    """
    Calculates the optimal poker strategy for the given scenario.
    
    The problem is to find the optimal strategy (action, sizing, frequency)
    for holding AA on a 2s2d2c board against a villain with KK.
    
    This is a GTO problem. With AA, we have the nuts. We need to decide
    whether to check or bet, and if betting, how much.

    1. The EV of betting is maximized by choosing the largest possible bet size.
       The EV function for a value bet of size B is: EV(B) = (B * Pot / (Pot + B)) + Pot
       With Pot=$10, EV(B) = 10*B/(10+B) + 10.
       This function increases as B increases, so the optimal bet size is the
       maximum possible, which is the stack size, $1000.

    2. We compare the EV of betting all-in vs. checking.
       EV(Bet 1000) = (1000 * 10 / (10 + 1000)) + 10 = 9.90 + 10 = 19.90
       EV(Check), where the villain may make a small value bet, is just over $10.
       Since EV(Bet) > EV(Check), we should always bet.

    3. The optimal strategy is therefore a pure strategy: bet 100% of the time.
    """
    
    # Define variables from the problem
    pot = 10
    stack = 1000
    
    # The optimal bet sizing is the largest possible bet, as it maximizes the EV of our value hand.
    bet_sizing = stack
    
    # Since betting is always higher EV than checking for our value hand, we do it 100% of the time.
    bet_frequency_percentage = 100
    
    # The problem asks to round sizing and percentage to the nearest even number.
    # 1000 is already an even number.
    final_sizing = 2 * round(bet_sizing / 2)
    
    # 100 is already an even number.
    final_frequency = 2 * round(bet_frequency_percentage / 2)

    # Format the final answer string.
    # The action is BET.
    action = "BET"
    
    # The "equation" is the final strategy string. We print the components.
    print(f"{action} {final_sizing} {final_frequency}%")

solve_poker_strategy()
# The final string is the answer to be parsed.
# The calculation shows a single optimal action.
# Action: BET
# Sizing: 1000
# Frequency: 100%
final_answer_string = "BET 1000 100%"
# The user asked me to return the answer with this specific format.
# Let's encapsulate the final string into the required format.
print(f"<<<{final_answer_string}>>>")