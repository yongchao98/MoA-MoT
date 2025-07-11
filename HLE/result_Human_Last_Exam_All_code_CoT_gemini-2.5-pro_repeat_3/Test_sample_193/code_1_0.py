import math

def solve_poker_strategy():
    """
    Calculates the optimal GTO strategy for Hero holding AA in the given scenario.
    """
    pot = 10
    stack = 1000
    
    # In a polarized vs. bluff-catcher scenario, we can model the equilibrium.
    # Hero's goal is to choose a bet size 'B' that maximizes his expected value (EV).
    # Villain's calling frequency 'q' will make Hero indifferent to bluffing: q = Pot / (Pot + B)
    # The EV for Hero's nut hand (AA) when betting B is: EV(AA) = Pot + q * B
    
    print("Step 1: Define the equation for the Expected Value (EV) of betting with AA.")
    print(f"The pot is ${pot}. The stack (maximum bet size) is ${stack}.")
    print("The EV equation for betting an amount 'B' with the nut hand (AA) is:")
    print("EV(AA) = Pot + (Pot / (Pot + B)) * B")
    print(f"EV(AA) = {pot} + ({pot} * B) / ({pot} + B)\n")

    # To maximize this EV, Hero needs to maximize B. The derivative of this EV function with
    # respect to B is always positive, meaning the function is always increasing with B.
    # Therefore, the optimal bet size is the largest possible bet.
    optimal_bet_size = stack
    
    print("Step 2: Find the optimal bet size 'B' that maximizes the EV.")
    print("The EV function increases as the bet size 'B' increases.")
    print(f"Therefore, the optimal bet size is the maximum possible amount, which is the stack: ${optimal_bet_size}\n")
    
    # The strategy for AA is to always take this highest EV action.
    action = "BET"
    sizing = optimal_bet_size
    frequency = 100
    
    print("Step 3: Determine the strategy for the AA hand.")
    print(f"The optimal action with AA is to bet ${sizing} with a frequency of {frequency}%.\n")

    # Rounding as per instructions: "Round both percentages and sizings to nearest even number"
    # 1000 is already an even number.
    if sizing % 2 != 0:
        sizing = int(2 * round(sizing / 2.0))
        
    # 100 is already an even number.
    if frequency % 2 != 0:
        frequency = int(2 * round(frequency / 2.0))
        
    print("Final Answer Formatting:")
    print(f"Action: {action}")
    print(f"Sizing: {sizing} (rounded to nearest even number)")
    print(f"Frequency: {frequency}% (rounded to nearest even number)\n")
    
    # The final answer in the required format
    final_answer_string = f"<<<{action} {sizing} {frequency}%>>>"
    
    # This is not part of the final answer, just showing the final equation components
    print("The final optimal strategy for AA is derived from this logic.")
    print(f"Final Action: {action}")
    print(f"Final Sizing: {sizing}")
    print(f"Final Frequency: {frequency}")


solve_poker_strategy()
