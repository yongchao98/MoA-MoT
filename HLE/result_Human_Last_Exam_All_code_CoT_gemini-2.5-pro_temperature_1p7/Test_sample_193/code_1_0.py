import math

def solve_poker_strategy():
    """
    This function calculates the optimal GTO strategy for holding AA
    in the described poker scenario. It explains the reasoning step-by-step.
    """
    
    # Step 1: Define the problem parameters based on the user's scenario.
    pot_size = 10
    stack_size = 1000
    action_type = "BET"

    print("--- Step 1: Problem Analysis ---")
    print(f"The pot is ${pot_size} and our effective stack is ${stack_size}.")
    print("Our hand is AA (the nuts). Our range is polarized (AA or QQ-33).")
    print("Villain's hand is KK, which acts as a bluff-catcher (it beats our bluffs but loses to our value hands).")
    print("This is a one-street betting game since future cards are blank.\n")

    print("--- Step 2: Finding the Optimal Bet Size ---")
    print("As the player out of position with a polarized range, our goal is to bet.")
    print("Our bet must be sized to maximize the value of our nut hand (AA).")
    print("The value from our AA bet comes when Villain calls. Villain will only call if we are also bluffing a certain percentage of the time.")
    print("Let 'B' be our bet size.")
    print("To make our bluffs indifferent between betting and checking (EV=0), Villain must call with a frequency 'alpha' where:")
    print("alpha = pot_size / (pot_size + B)")
    print("The profit we make from our AA bet is `profit = alpha * B`.")
    print(f"Substituting 'alpha', we want to maximize the function: profit(B) = (pot_size * B) / (pot_size + B)")
    print(f"For our scenario, this is profit(B) = ({pot_size} * B) / ({pot_size} + B).")
    print("This function's value increases as 'B' increases. To maximize our profit, we must choose the largest possible bet size.\n")
    
    optimal_sizing = stack_size
    print(f"The largest possible bet is our entire stack. Therefore, the optimal bet sizing is: ${optimal_sizing}\n")

    print("--- Step 3: Finding the Optimal Betting Frequency for AA ---")
    print("We have determined that betting the maximum size ($1000) maximizes the EV of our value hands.")
    print("If we check with AA, we risk the villain checking behind with KK, and we would only win the small starting pot.")
    print("Betting our entire stack forces the villain into a tough spot and extracts the maximum possible value when he calls.")
    print("Therefore, the optimal strategy is to bet with AA 100% of the time.\n")
    
    bet_frequency = 100

    # Step 4: Round values to the nearest even number and format the output.
    final_sizing = 2 * round(optimal_sizing / 2)
    final_frequency = 2 * round(bet_frequency / 2)

    print("--- Step 4: Final Optimal Strategy ---")
    print(f"The action is {action_type}.")
    print(f"The bet sizing is {final_sizing}.")
    print(f"The frequency is {final_frequency}%.")

    print("\nFinal Answer:")
    # We must output each number that forms the final answer string.
    print(f"{action_type} {final_sizing} {final_frequency}%")


solve_poker_strategy()