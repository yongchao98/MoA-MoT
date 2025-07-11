import math

def find_optimal_strategy():
    """
    Solves for the optimal poker strategy in the given scenario.
    """
    # Step 1: Establish the Game State
    pot = 10
    stack = 1000

    print("Step-by-step analysis of the optimal strategy for Hero with AA:")
    print(f"Pot (P) = ${pot}, Effective Stack = ${stack}\n")

    # Step 2 & 3: Analyze the EV of Checking with AA
    print("--- Analyzing Option 1: CHECK with AA ---")
    print("If Hero checks, Villain (holding KK) must decide whether to bet or check.")
    print("Villain knows Hero's range consists of AA and QQ-33. A bet from Villain would only be called or raised by AA; all of Hero's other hands would fold.")
    print("Since a bet only gets action from a better hand, Villain's optimal play is to check behind.")
    print("When Villain checks behind, the hand goes to showdown, and Hero's AA wins the pot.")
    
    ev_check = pot
    print("\nThe Expected Value (EV) of checking with AA is the value of the pot.")
    print(f"EV(Check AA) = P = ${ev_check}\n")

    # Step 4: Analyze the EV of Betting with AA
    print("--- Analyzing Option 2: BET 'B' with AA ---")
    print("For a bet with AA to get value, Hero must balance it by bluffing with some weaker hands (QQ-33).")
    print("A GTO strategy requires that Hero's bluffing frequency makes Villain indifferent to calling or folding with his KK.")
    print("This, in turn, dictates Villain's optimal calling frequency 'c'. To make Hero indifferent to bluffing (i.e., EV of bluffing = 0), Villain must call with a frequency of:")
    print("c = Pot / (Pot + Bet)\n")
    
    print("Given Villain's calling frequency 'c', we can calculate the EV of betting with AA.")
    print("EV(Bet AA) = (c * Value_if_Called) + ((1-c) * Value_if_Folds)")
    print("EV(Bet AA) = c * (Pot + Bet) + (1-c) * Pot")
    print("By substituting 'c', this simplifies to:")
    print("EV(Bet AA) = Pot + (Bet / (Pot + Bet)) * Pot\n")

    # Step 5: Compare EVs and Determine the Optimal Bet Size
    print("--- Comparing Strategies and Finding the Optimal Bet Size ---")
    print("We compare the EV of our two options for a given bet size 'B':")
    print(f"EV(Bet AA) = {pot} + ({pot} * B) / ({pot} + B)")
    print(f"EV(Check AA) = {pot}\n")
    
    print("Since the term '(B / (Pot + B)) * Pot' is always positive for any B > 0, the EV of Betting is always greater than the EV of Checking.")
    print("Therefore, Hero should always bet with AA. The frequency is 100%.")

    print("\nTo maximize profit, Hero must choose the bet size 'B' that maximizes the EV(Bet AA) function.")
    print("The value of the function increases as 'B' increases. Thus, the optimal bet size is the largest possible bet: the entire stack.")
    
    optimal_bet_size = stack
    ev_bet_optimal = pot + (pot * optimal_bet_size) / (pot + optimal_bet_size)

    print(f"\nFinal EV calculation with the optimal bet size B = ${optimal_bet_size}:")
    print(f"EV(Bet AA) = {pot} + ({optimal_bet_size} / ({pot} + {optimal_bet_size})) * {pot} = ${ev_bet_optimal:.2f}\n")

    # Step 6: Format the Final Answer
    action = "BET"
    sizing = optimal_bet_size
    frequency = 100

    # Round to the nearest even number per instructions
    final_sizing = int(round(sizing / 2.0)) * 2
    final_frequency = int(round(frequency / 2.0)) * 2
    
    final_strategy_string = f"{action} {final_sizing} {final_frequency}%"

    print("--- Conclusion: Optimal Strategy for AA ---")
    print(f"The final strategy is: {final_strategy_string}")

find_optimal_strategy()